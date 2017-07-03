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
	private HashMap<String, Long> contigPointers;
	
	private int chromCol = -1;
	private int startCol = -1;
	private int refCol = -1;
	private int altCol = -1;

	Set<Integer> startSet = new HashSet<Integer>();
	ArrayList<VariantContext> possibleVariants = new ArrayList<VariantContext>();

	public FilteredVariantReader(File inFile) throws Exception {
		try {
			raFile = new RandomAccessFile(inFile, "r");
			contigPointers = new HashMap<String, Long>();

			String line = null;
			
			// TODO: If file ends with .gz, handle!!!
			if(inFile.getPath().endsWith(".vcf") || inFile.getPath().endsWith(".vcf.gz")){				
				// Skip all commented lines if VCF
				while ((line = raFile.readLine()).startsWith("##")) {
				}
				
				if(!line.startsWith("#")){
					throw new Error("Invalid vcf. Could not find header.");
				} else {
					line = line.substring(1);
				}
				
				// Vcf is standard
				chromCol = 0;
				startCol = 1;
				refCol = 3;
				altCol = 4;
			} else {
				// If not vcf, assume gemini and attempt to parse
				line = raFile.readLine();
				
				// Parse col locations
				String[] headers = line.split("\t");
				for(int x = 0; x < headers.length; x++){
					if(headers[x].equals("chrom")){
						chromCol = x;
					} else if(headers[x].equals("start")){
						startCol = x;
					} else if(headers[x].equals("ref")){
						refCol = x;
					} else if(headers[x].equals("alt")){
						altCol = x;
					}
				}
				if(chromCol == -1 || startCol == -1 || refCol == -1 || altCol == -1){
					throw new Exception("Filtered variants weren't provided in either vcf or gemini format. No other format is supported. Please supply a vcf or gemini output file.");
				}
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
				String line = raFile.readLine();
				
				if(line != null){
					String[] entries = line.split("\\t");

					// Parse all alleles
					ArrayList<Allele> alleles = new ArrayList<Allele>();
					Allele allele = Allele.create(entries[refCol], true);
					alleles.add(allele);
					String[] nonRefAlleles = entries[altCol].split(",");
					for (String alleleString : nonRefAlleles) {
						Allele a = Allele.create(alleleString, false);
						alleles.add(a);
					}

					long stop;
					long start = Long.parseLong(entries[startCol]);
					stop = start + allele.length() - 1;
					

					// Create new variantContext and add
					if (!startSet.contains(entries.hashCode())) {
						possibleVariants.add(
								vcBuilder.source(fileName).chr(entries[chromCol]).start(start).stop(stop).alleles(alleles).make());
					}
					
					startSet.add(entries.hashCode());
					
					// Check if contig changed
					if (!possibleVariants.get(possibleVariants.size() - 1).getContig().equals(curInterval.getContig())) {
						break;
					}
				} else {
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
