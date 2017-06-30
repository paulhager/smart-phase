import java.io.EOFException;
import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import htsjdk.samtools.util.Interval;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class FilteredVariantReader {

	private RandomAccessFile raFile;
	private String fileName;
	private boolean vcf = false;
	private HashMap<String, Long> contigPointers;

	Set<Integer> startSet = new HashSet<Integer>();
	ArrayList<VariantContext> possibleVariants = new ArrayList<VariantContext>();

	public FilteredVariantReader(File inFile) {
		try {
			raFile = new RandomAccessFile(inFile, "r");
			contigPointers = new HashMap<String, Long>();

			String line;

			// Skip all commented lines if VCF
			
			while ((line = raFile.readLine()).startsWith("##")) {
			}

			// Remove comment char from VCF header
			if (line.startsWith("#")) {
				vcf = true;
				line = line.substring(1);
			}

			// Ensure correct header
			String[] header = line.split("\\t");
			if (!header[0].equalsIgnoreCase("chrom")
					|| !(header[1].equalsIgnoreCase("start") || header[1].equalsIgnoreCase("pos"))) {
				throw new Error("Incorrect header in filtered variant file.");
			}
			
			// Save all contig start positions in file
			String curContig = "";
			long prevLinePointer = raFile.getFilePointer();
			while(true) {
				try {
					line = raFile.readLine();
					if(line == null) {
						break;
					}
					if(!line.split("\\t")[0].equals(curContig)) {
						curContig = line.split("\\t")[0];
						contigPointers.put(curContig, prevLinePointer);
					}
					prevLinePointer = raFile.getFilePointer();
				} catch (EOFException e) {
					break;
				}
			}
		} catch (Error | IOException e) {
			e.printStackTrace();
		}
	}

	public void close() {
		try {
			raFile.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	public boolean contigImportantCheck(String contig) {
		return contigPointers.containsKey(contig);
	}

	public ArrayList<VariantContext> scan(Interval curInterval, boolean contigSwitch) throws IOException {

		int intervalStart = curInterval.getStart();
		int intervalEnd = curInterval.getEnd();
		
		// Contig changed, so delete all variants not on new contig 
		if(contigSwitch){
			possibleVariants = new ArrayList<VariantContext>();
			raFile.seek(contigPointers.get(curInterval.getContig()));
		}

		VariantContextBuilder vcBuilder = new VariantContextBuilder();
		
		// Add variants until variant exceeds intervalEnd or reader is empty
		while ((possibleVariants.size() == 0
				|| possibleVariants.get(possibleVariants.size() - 1).getStart() < intervalEnd)) {
			try {
				String[] entries = raFile.readLine().split("\\t");

				// Parse all alleles
				ArrayList<Allele> alleles = new ArrayList<Allele>();
				Allele allele = Allele.create(entries[3], true);
				alleles.add(allele);
				String[] nonRefAlleles = entries[4].split(",");
				for (String alleleString : nonRefAlleles) {
					Allele a = Allele.create(alleleString, false);
					alleles.add(a);
				}

				long stop;
				long start = Long.parseLong(entries[1]);
				if (vcf) {
					stop = start + allele.length() - 1;
				} else {
					stop = Long.parseLong(entries[2]);
				}

				// Create new variantContext and add
				if (!startSet.contains(entries.hashCode())) {
					possibleVariants.add(
							vcBuilder.source(fileName).chr(entries[0]).start(start).stop(stop).alleles(alleles).make());
				}
				
				startSet.add(entries.hashCode());
				
				// Check if contig changed
				if (!possibleVariants.get(possibleVariants.size() - 1).getContig().equals(curInterval.getContig())) {
					break;
				}
			} catch (EOFException e) {
				break;
			}
			
		}

		// Remove all variants less than current interval start
		possibleVariants.removeIf(v -> v.getStart() < intervalStart && v.getContig().equals(curInterval.getContig()));

		// Create to-be-returned variants arraylist
		// TODO: PERFORMANCE BOOST IF YOU CHECK FROM BEHIND?
		ArrayList<VariantContext> variants = new ArrayList<>(possibleVariants);
		variants.removeIf(v -> v.getStart() > intervalEnd);

		return variants;
	}

}
