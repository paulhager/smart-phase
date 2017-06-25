import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

import htsjdk.samtools.util.Interval;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class FilteredVariantReader {
	
	private BufferedReader fileReader;
	private String fileName;
	private boolean vcf = false;
	Set<Long> startSet = new HashSet<Long>();
	ArrayList<VariantContext> possibleVariants = new ArrayList<VariantContext>();
	
	public FilteredVariantReader(File inFile){
		try {
			fileName = inFile.getName();
			fileReader = new BufferedReader(new FileReader(inFile));
			
			String line;
			
			// Skip all commented lines if VCF
			while((line = fileReader.readLine()).startsWith("##")){}
			
			// Remove comment char from VCF header
			if(line.startsWith("#")){
				vcf = true;
				line = line.substring(1);
			}

			// Ensure correct header 
			String[] header = line.split("\\t");
			if(!header[0].equalsIgnoreCase("chrom") || !(header[1].equalsIgnoreCase("start") || header[1].equalsIgnoreCase("pos"))){
				throw new Error("Incorrect header in filtered variant file.");
			}
		} catch (Error | IOException e) {
			e.printStackTrace();
		}
	}
	
	public void close(){
		try {
			fileReader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public ArrayList<VariantContext> scan(Interval curInterval) throws IOException{
		
		int intervalStart = curInterval.getStart();
		int intervalEnd = curInterval.getEnd();
		
		VariantContextBuilder vcBuilder = new VariantContextBuilder(); 
		
		// Add variants until variant exceeds intervalEnd or reader is empty
		while((possibleVariants.size() == 0 || possibleVariants.get(possibleVariants.size()-1).getStart() < intervalEnd) && fileReader.ready()){
			String[] entries = fileReader.readLine().split("\\t");
			
			// Parse all alleles
			ArrayList<Allele> alleles = new ArrayList<Allele>();
			Allele allele = Allele.create(entries[3], true);
			alleles.add(allele);
			String[] nonRefAlleles = entries[4].split(",");
			for(String alleleString : nonRefAlleles){
				Allele a = Allele.create(alleleString, false);
				alleles.add(a);
			}
			
			long stop;
			long start = Long.parseLong(entries[1]);
			if(vcf){
				stop = start + allele.length() - 1;
			} else {
				stop = Long.parseLong(entries[2]);
			}
			
			// Create new variantContext and add
			if(!startSet.contains(start)){
				possibleVariants.add(vcBuilder.source(fileName).chr(entries[0]).start(start).stop(stop).alleles(alleles).make());
			}
			startSet.add(start);
		}
		
		// Remove all variants less than current interval start
		possibleVariants.removeIf(v -> v.getStart() < intervalStart);
		
		
		// Create to-be-returned variants arraylist
		// TODO: PERFORMANCE BOOST IF YOU CHECK FROM BEHIND?
		ArrayList<VariantContext> variants = new ArrayList<>(possibleVariants);
		variants.removeIf(v -> v.getStart() > intervalEnd);
		
		return variants;
	}

}
