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
	
	private String spliter = "\\t";
	
	private boolean vcf = false;
	private boolean cohort;
	
	private ArrayList<VariantContext> allVariants = new ArrayList<VariantContext>();

	Set<String[]> startSet = new HashSet<String[]>();
	Set<String> allVarsContigs = new HashSet<String>();
	ArrayList<VariantContext> possibleVariants = new ArrayList<VariantContext>();

	public FilteredVariantReader(File inFile, boolean cohort, String patID) throws Exception {
		this.cohort = cohort;
		try {
			raFile = new RandomAccessFile(inFile, "r");
			contigPointers = new HashMap<String, Long>();

			String line = null;
			
			if(inFile.getPath().endsWith(".tsv")){
				spliter = ",";
			}
			
			if(this.cohort){
				line = raFile.readLine();
				while(true) {
					try {
						line = raFile.readLine();
						if(line == null) {
							break;
						}
						String[] cols=line.split(",");
						if(cols[0].indexOf(patID+" ") != -1){
							String var1 = cols[1];
							String var2 = cols[2];
							
							String[] var1Data = var1.split("-");
							String[] var2Data = var2.split("-");
							
							//Ugly but im lazy
							
							VariantContext var1VC = createVC(var1Data);
							VariantContext var2VC = createVC(var2Data);
							
							allVariants.add(var1VC);
							allVariants.add(var2VC);
							
							allVarsContigs.add(var1VC.getContig());
							
						}
					} catch (EOFException e) {
						return;
					}
				}
				return;
			}
			
			// TODO: If file ends with .gz, handle!!!
			if(inFile.getPath().endsWith(".vcf") || inFile.getPath().endsWith(".vcf.gz")){		
				vcf = true;
				
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
				String[] headers = line.split(spliter);
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
					if(!line.split(spliter)[0].equals(curContig)) {
						curContig = line.split(spliter)[0];
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
	
	private VariantContext createVC(String[] data){
		// Parse all alleles
		ArrayList<Allele> alleles = new ArrayList<Allele>();
		Allele allele = Allele.create(data[2], true);
		alleles.add(allele);
		String[] nonRefAlleles = data[3].split(",");
		for (String alleleString : nonRefAlleles) {
			Allele a = Allele.create(alleleString, false);
			alleles.add(a);
		}

		long stop;
		long start = Long.parseLong(data[1]);
		
		stop = start + allele.length() - 1;
		
		return new VariantContextBuilder().source(fileName).chr(data[0]).start(start).stop(stop).alleles(alleles).make();
	}
	
	
	public boolean contigImportantCheck(String contig) {
		if(cohort){
			return allVarsContigs.contains(contig);
		}
		
		return contigPointers.containsKey(contig);
	}
	
	public ArrayList<VariantContext> scan(Interval curInterval, boolean contigSwitch) throws IOException {

		int intervalStart = curInterval.getStart();
		int intervalEnd = curInterval.getEnd();
		
		if(cohort){
			ArrayList<VariantContext> scanVars = new ArrayList<VariantContext>(allVariants);
			scanVars.removeIf(v -> v.getContig() != curInterval.getContig() || v.getStart() < intervalStart || v.getEnd() > intervalEnd);
			return scanVars;
		}
		
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
					String[] entries = line.split(spliter);

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
					
					// Gemini is 0-based but rest of program assumes 1-based
					if(!vcf){
						start++;
					}
					stop = start + allele.length() - 1;

					// Create new variantContext and add
					if (!startSet.contains(entries)) {
						possibleVariants.add(
								vcBuilder.source(fileName).chr(entries[chromCol]).start(start).stop(stop).alleles(alleles).make());
					}
					
					startSet.add(entries);
					
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
