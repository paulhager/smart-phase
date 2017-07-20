import java.io.EOFException;
import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFFileReader;

public class FilteredVariantReader {

	private RandomAccessFile raFile;
	private String fileName;
	private HashMap<String, Long> contigPointers;
	
	private int chromCol = -1;
	private int startCol = -1;
	private int refCol = -1;
	private int altCol = -1;
	
	private String spliter = "\\t";
	
	private String patientID = "";
	
	private boolean vcf = false;
	private boolean cohort;
	private boolean gzipVCF = false;
	
	private VCFFileReader vcfREADER;
	
	private ArrayList<VariantContext> allVariants = new ArrayList<VariantContext>();

	Set<String> allVarsContigs = new HashSet<String>();
	ArrayList<VariantContext> possibleVariants = new ArrayList<VariantContext>();
	
	IntervalList iList;

	public FilteredVariantReader(File inFile, boolean cohort, String patID, IntervalList iList) throws Exception {
		patientID = patID;
		fileName = inFile.getName();
		this.iList = iList;
		this.cohort = cohort;
		try {
			raFile = new RandomAccessFile(inFile, "r");
			contigPointers = new HashMap<String, Long>();

			String line = null;
			
			if(inFile.getPath().endsWith(".csv")){
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
						String[] cols=line.split(spliter);
						if(cols.length==1){
							throw new Exception("While parsing filtered variants, tried to split using \'"+spliter+"\' but failed. Please ensure your file is tab-seperated, or ends in .csv if it is comma separated.");
						}
						if(cols.length == 2){
							continue;
						}
						List<String> samples = Arrays.asList(cols[0].split(","));
						for(String sample : samples){
							samples.set(samples.indexOf(sample), sample.trim());
						}
						if(samples.contains(patID)){
							String var1 = cols[1];
							String var2 = cols[2];
							
							String[] var1Data = var1.split("-");
							String[] var2Data = var2.split("-");
														
							createVC(var1Data);
							createVC(var2Data);
							
							String chrom = var1Data[0];
							if(!chrom.startsWith("chr")){
								chrom = "chr"+chrom;
							}
							allVarsContigs.add(chrom);
							
							int intStart = Integer.parseInt(var1Data[1]);
							intStart = intStart - 100;
							int intEnd = Integer.parseInt(var2Data[1]);
							intEnd = intEnd + 100;
							
							this.iList.add(new Interval(chrom, intStart, intEnd));
						}
					} catch (EOFException e) {
						return;
					}
				}
				return;
			}
			
			//File is bgzipped vcf file and can be handeled as such. Index is required
			if(inFile.getPath().endsWith(".vcf.gz")){
				vcfREADER = new VCFFileReader(inFile);
				gzipVCF = true;
				return;
			}
			
			// File is normal vcf
			if(inFile.getPath().endsWith(".vcf")){		
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
					if(headers[x].trim().equals("chrom")){
						chromCol = x;
					} else if(headers[x].trim().equals("start")){
						startCol = x;
					} else if(headers[x].trim().equals("ref")){
						refCol = x;
					} else if(headers[x].trim().equals("alt")){
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
			int prevStart = 0;
			while(true) {
				try {
					line = raFile.readLine();
					if(line == null) {
						break;
					}
					
					String[] splitLine = line.split(spliter);
					String actContig = splitLine[chromCol];
					if(!actContig.equals(curContig)) {
						curContig = actContig;
						if(contigPointers.containsKey(curContig)){
							throw new Exception("ERROR! Filtered variants list must be grouped by chromosome and sorted by start position within each chromosome.");
						}
						contigPointers.put(curContig, prevLinePointer);
						prevStart = 0;
					}
					int curStart = Integer.parseInt(splitLine[startCol]);
					// Ensure within chr. that all vars are start sorted
					if(prevStart > curStart){
						throw new Exception("ERROR! Filtered variants list must be grouped by chromosome and sorted by start position within each chromosome. The offending variants are: (prev) "+actContig+"-"+prevStart+" || "+actContig+"-"+curStart);
					}
					prevStart = curStart;
					
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
	
	public IntervalList getiList(){
		return this.iList;
	}
	
	private void createVC(String[] data){
		// Parse all alleles
		ArrayList<Allele> alleles = new ArrayList<Allele>();
		Allele allele = Allele.create(data[2], true);
		alleles.add(allele);
		Allele altAllele = Allele.create(data[3], false);
		alleles.add(altAllele);

		long stop;
		long start = Long.parseLong(data[1]);
		
		stop = start + allele.length() - 1;
		
		String chrom = data[0];
		if(!chrom.startsWith("chr")){
			chrom = "chr"+chrom;
		}
		
		
		
		VariantContext varVC = new VariantContextBuilder().source(fileName).chr(chrom).start(start).stop(stop).alleles(alleles).genotypes(new GenotypeBuilder().alleles(alleles).name(patientID).make()).make();
		
		boolean notNew = false;
		// Check if variant already exists
		for(VariantContext posVar : allVariants){
			if(varVC.getStart() == posVar.getStart() && varVC.getEnd() == posVar.getEnd() && varVC.hasSameAllelesAs(posVar) && varVC.hasSameAlternateAllelesAs(posVar) && varVC.getContig().equals(posVar.getContig())){
				notNew = true;
			}
		}
		
		if(!notNew){
			allVariants.add(varVC);
		} 
		
		return;
	}
	
	
	public boolean contigImportantCheck(String contig) {
		if(cohort){
			return allVarsContigs.contains(contig);
		}
		
		if(gzipVCF){
			return true;
		}
		
		return contigPointers.containsKey(contig);
	}
	
	public ArrayList<VariantContext> scan(Interval curInterval, boolean contigSwitch) throws IOException {

		int intervalStart = curInterval.getStart();
		int intervalEnd = curInterval.getEnd();
		
		if(gzipVCF){
			return new ArrayList<VariantContext>(vcfREADER.query(curInterval.getContig(), intervalStart, intervalEnd).toList());
		}
		
		if(cohort){
			System.out.println(allVariants.size());
			ArrayList<VariantContext> scanVars = new ArrayList<VariantContext>(allVariants);
			scanVars.removeIf(v -> !v.getContig().equals(curInterval.getContig()) || v.getStart() < intervalStart || v.getEnd() > intervalEnd);
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
				|| (possibleVariants.get(possibleVariants.size() - 1).getStart() < intervalEnd && possibleVariants.get(possibleVariants.size() - 1).getContig().equals(curInterval.getContig())))) {
			try {
				String line = raFile.readLine();
				
				if(line != null){
					String[] entries = line.split(spliter);

					// Parse all alleles
					ArrayList<Allele> alleles = new ArrayList<Allele>();
					Allele allele = Allele.create(entries[refCol], true);
					alleles.add(allele);
					
					// Ensure file is normalized
					String nonRefAllele = entries[altCol];
					if(nonRefAllele.indexOf(",") != -1){
						throw new Exception("Only normalized variant containing files are accepted!");
					}
					Allele a = Allele.create(nonRefAllele, false);
					alleles.add(a);
					
					

					long stop;
					long start = Long.parseLong(entries[startCol]);
					
					// Gemini is 0-based but rest of program assumes 1-based
					if(!vcf){
						start++;
					}
					stop = start + allele.length() - 1;
					
					VariantContext newVarC = vcBuilder.source(fileName).chr(entries[chromCol]).start(start).stop(stop).alleles(alleles).genotypes(new GenotypeBuilder().alleles(alleles).name(patientID).make()).make();

					boolean notNew = false;
					// Check if variant already exists
					for(VariantContext posVar : possibleVariants){
						if(newVarC.getStart() == posVar.getStart() && newVarC.getEnd() == posVar.getEnd() && newVarC.hasSameAllelesAs(posVar) && newVarC.hasSameAlternateAllelesAs(posVar) && newVarC.getContig().equals(posVar.getContig())){
							notNew = true;
						}
					}
					// Create new variantContext and add
					if (!notNew && curInterval.getContig().equals(entries[chromCol])) {
						possibleVariants.add(newVarC);
					} 					
					
					
					// Check if contig changed
					if (!newVarC.getContig().equals(curInterval.getContig())) {
						break;
					}
				} else {
					break;
				}
				
			} catch (EOFException e) {
				break;
			} catch (Exception e) {
				e.printStackTrace();
				System.exit(1);
			}
			
		}

		// Remove all variants less than current interval start
		possibleVariants.removeIf(v -> v.getStart() < intervalStart || !v.getContig().equals(curInterval.getContig()));

		// Create to-be-returned variants arraylist
		// TODO: PERFORMANCE BOOST IF YOU CHECK FROM BEHIND?
		ArrayList<VariantContext> variants = new ArrayList<>(possibleVariants);
		variants.removeIf(v -> v.getStart() > intervalEnd);

		return variants;
	}

}
