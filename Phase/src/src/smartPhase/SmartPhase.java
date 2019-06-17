//  Copyright (C) 2018 the SmartPhase contributors.
//  Website: https://github.com/paulhager/smart-phase
//
//  This file is part of the SmartPhase phasing tool.
//
//  The SmartPhase phasing tool is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

package smartPhase;

import java.io.BufferedReader;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.ListIterator;
import java.util.Set;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;

//import org.apache.commons.lang.*;
//import org.broadinstitute.gatk.utils.GenomeLocSortedSet.MergeStrategy;
//import org.broadinstitute.gatk.utils.codecs.samread.SAMReadCodec;
//import org.broadinstitute.gatk.utils.genotyper.AlleleList;
import org.apache.commons.cli.*;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import picard.pedigree.PedFile;
import smartPhase.HaplotypeBlock.Strand;

public class SmartPhase {

	public enum Phase {
		CIS, TRANS, TOTAL_OBSERVED
	}

	// static SAMRecordIterator samIterator;
	static ArrayList<SAMRecord> curRecords = new ArrayList<SAMRecord>();
	static ArrayList<SAMRecordIterator> samIteratorList = new ArrayList<SAMRecordIterator>();
	static HashMap<SAMRecordIterator, SAMRecord> grabLastRec = new HashMap<SAMRecordIterator, SAMRecord>();
	static HashMap<String, SAMRecord> pairedEndReads = new HashMap<String, SAMRecord>();
	static HashMap<String, SAMRecord> pairedEndReadsHelperRecord = new HashMap<String, SAMRecord>();
	static HashMap<String, SAMRecord> pairedEndReadsHelperName = new HashMap<String, SAMRecord>();
	static SAMRecord curRec = null;
	static File[] inputREADFILES = null;
	static String prevContig = "";
	static double[] minMAPQ;
	static double vcfCutoff = 0;

	static HashSet<VariantContext> seenInRead = new HashSet<VariantContext>();
	static HashSet<VariantContext> NOT_SeenInRead = new HashSet<VariantContext>();
	

	static HashMap<VariantContext, Integer> varToStartHash = new HashMap<VariantContext, Integer>();
	static HashMap<VariantContext, Integer> varToEndHash = new HashMap<VariantContext, Integer>();

	static HashSet<VariantContext> neverSeenVariants = new HashSet<VariantContext>();

	static HashMap<PhaseCountTriple<Set<VariantContext>, Phase>, Double> phaseCounter = new HashMap<PhaseCountTriple<Set<VariantContext>, Phase>, Double>();
	static HashMap<VariantContext, VariantContext> pairedReadsVarMaps = new HashMap<VariantContext, VariantContext>();

	static String PATIENT_ID;

	static boolean TRIO = false;
	static boolean READS = false;
	static boolean PHYSICAL_PHASING = false;
	static boolean REJECT_PHASE = false;
	static boolean VALIDATION = false;
	static boolean WRITE_VCF = false;

	static boolean countReads = false;

	// Statistics
	static int denovoCounter = 0;
	static int notSimpleCounter = 0;
	static int totalVarCounter = 0;
	static int globalCis = 0;
	static int globalCisLength = 0;
	static int globalTrans = 0;
	static int globalTransLength = 0;
	static int globalNewBlock = 0;
	static int globalContradiction = 0;
	static int globalRNAseqCount = 0;
	static int ppPhased = 0;

	static HashMap<String, String> physicalPhasingPIDMap = new HashMap<String, String>();
	static HashMap<String, String> physicalPhasingPGTMap = new HashMap<String, String>();

	static HashMap<PhaseCountTriple<Set<VariantContext>, Phase>, Double> skipIntronCounter = new HashMap<PhaseCountTriple<Set<VariantContext>, Phase>, Double>();

	static HashMap<VariantContext, VariantContext> exonStartVars = new HashMap<VariantContext, VariantContext>();

	public static void main(String[] args) throws Exception {
		// Parse Options
		ArrayList<Option> optionsList = new ArrayList<>();
		
		Option filteredVariantsOption = Option.builder("f").longOpt("filtered-variants").argName("file").hasArg().desc("Path to file containing patient variants filtered for significance").build();
		optionsList.add(filteredVariantsOption);
		
		Option allVariantsOption = Option.builder("a").longOpt("all-variants").argName("file.vcf").hasArg().required().desc("Path to file containing all patient variants").build();
		optionsList.add(allVariantsOption);
		
		Option patientOption = Option.builder("p").longOpt("patient").argName("string").hasArg().required().desc("ID of patient used in vcf and ped files").build();
		optionsList.add(patientOption);
		
		Option outputOption = Option.builder("o").longOpt("output").argName("file").hasArg().required().desc("Path to desired output file").build();
		optionsList.add(outputOption);
		
		Option geneRegionsOption = Option.builder("g").longOpt("gene-regions").argName("file.bed").hasArg().desc("Path to file containing genomic regions to be analyzed").build();
		optionsList.add(geneRegionsOption);
		
		Option readsOption = Option.builder("r").longOpt("reads").argName("file.bam,file.bam,...").hasArgs().desc("Comma seperated list of paths to files containing aligned patient reads").build();
		optionsList.add(readsOption);
		
		Option mapqOption = Option.builder("m").longOpt("mapq").argName("int,int,...").hasArgs().desc("Comma seperated list of mapping quality cutoff values to use when examining reads. Each value corresponds to the min MAPQ for an input BAM file").build();
		optionsList.add(mapqOption);
				
		Option trioOption = Option.builder("t").longOpt("trio").desc("Indicates trio information is available and contained in original-variants file provided").build();
		optionsList.add(trioOption);
		
		Option pedOption = Option.builder("d").longOpt("ped").argName("file.ped").hasArg().desc("Path to file containing vcf IDs of trio").build();
		optionsList.add(pedOption);
		
		Option physicalPhasingOption = Option.builder("y").longOpt("physical-phasing").desc("Indicates GATK physical phasing info should be used as last-resort phasing").build();
		optionsList.add(physicalPhasingOption);
		
		Option rejectPhaseOption = Option.builder("x").longOpt("reject-phase").desc("Indicates phase already present in vcf should be disregarded").build();
		optionsList.add(rejectPhaseOption);
		
		Option validationOption = Option.builder("v").longOpt("validation").desc("Internal validation option. Do not use.").build();
		optionsList.add(validationOption);
		
		Option vcfWriteOption = Option.builder("vcf").desc("Indicates a vcf file should be generated to include the results of phasing. Requires a confidence cutoff specified through the --confidence option.").build();
		optionsList.add(vcfWriteOption);
		
		Option vcfWriteCutoff = Option.builder("c").longOpt("cutoff").hasArg().desc("The confidence cutoff to be used when determining if a variant should be indicated as phased in the vcf. Only used in conjunction with the --vcf option.").build();
		optionsList.add(vcfWriteCutoff);
		
		Option helpOption = Option.builder("h").longOpt("help").desc("Print this message").build();
		optionsList.add(helpOption);
				
		Options options = new Options();
		Options helpCheckOptions = new Options();
		
		for(Option option : optionsList) {
			options.addOption(option);
			boolean required = option.isRequired();
			option.setRequired(false);
			helpCheckOptions.addOption(option);
			option.setRequired(required);
		}
		
		CommandLineParser parser = new DefaultParser();
		CommandLine helpCheckCmd = parser.parse(helpCheckOptions, args);
					
		if(helpCheckCmd.hasOption(helpOption.getLongOpt()) || helpCheckCmd.hasOption(helpOption.getOpt())) {
			HelpFormatter helpFormatter = new HelpFormatter();
			String usageMessage = "Welcome to SmartPhase! A dedicated tool designed to assist in the rapid and accurate phasing of variant combinations for clinical analysis. Please refer to the following list of options on how to pass the necessary parameters for use:\n";
			helpFormatter.printHelp(usageMessage, options, true);
			System.exit(0);
		}
		CommandLine cmd = null;
		try {
			cmd = parser.parse(options, args);			
		} catch (ParseException pe) {
			HelpFormatter helpFormatter = new HelpFormatter();
			String usageMessage = "Welcome to SmartPhase! A dedicated tool designed to assist in the rapid and accurate phasing of variant combinations for clinical analysis. Please refer to the following list of options on how to pass the necessary parameters for use:\n";
			helpFormatter.printHelp(usageMessage, options, true);
			System.exit(1);
		}

		long startTime = System.currentTimeMillis();

		File inputBED = null;
		final File inputVCF_ALL = new File(cmd.getOptionValue(allVariantsOption.getOpt()));
		final String inputREADFILESSTRING = cmd.getOptionValue(readsOption.getOpt());
		final String inputMinMAPQ = cmd.getOptionValue(mapqOption.getOpt());
		final Path OUTPUT = Paths.get(cmd.getOptionValue(outputOption.getOpt()));
		PATIENT_ID = cmd.getOptionValue(patientOption.getOpt());
		final boolean PAIRED;
		final File inputPEDIGREE;
		PedFile familyPed = null;
		final File inputVCF_FILTER;
		Set<String> unimportantContigs = new HashSet<>();
		ArrayList<String> noPrintAttributes = new ArrayList<String>();
		noPrintAttributes.add("Preceding");
		noPrintAttributes.add("Innocuous");
		noPrintAttributes.add("mergedBlocks");
		noPrintAttributes.add("linkedPreceding");
		noPrintAttributes.add("linkedConfidence");
		
		if(cmd.hasOption(filteredVariantsOption.getOpt())) {
			inputVCF_FILTER = new File(cmd.getOptionValue(filteredVariantsOption.getOpt()));
		} else {
			inputVCF_FILTER = inputVCF_ALL;
		}
		
		if (cmd.hasOption(trioOption.getOpt())) {
			TRIO = true;
			inputPEDIGREE = new File(cmd.getOptionValue(pedOption.getOpt()));
			familyPed = PedFile.fromFile(inputPEDIGREE, true);
		}

		if (cmd.hasOption(rejectPhaseOption.getOpt())) {
			REJECT_PHASE = true;
		}

		if (cmd.hasOption(physicalPhasingOption.getOpt())) {
			PHYSICAL_PHASING = true;
		}

		if (cmd.hasOption(validationOption.getOpt())) {
			VALIDATION = true;
		}
		
		if(cmd.hasOption(vcfWriteOption.getOpt())) {
			if(!cmd.hasOption(vcfWriteCutoff.getOpt())) {
				throw new Exception("If a vcf file with the phasing results should be generated (-vcf), a confidence cutoff must be provided (-c)");				
			} else {
				WRITE_VCF = true;
				vcfCutoff = Float.parseFloat(cmd.getOptionValue(vcfWriteCutoff.getOpt()));				
			}
		}

		if (cmd.hasOption(geneRegionsOption.getOpt())) {
			PAIRED = false;

			inputBED = new File(cmd.getOptionValue(geneRegionsOption.getOpt()));

			if (!inputBED.exists()) {
				throw new FileNotFoundException("File " + inputBED.getAbsolutePath() + " does not exist!");
			}
		} else {
			PAIRED = true;
		}
	

		if (cmd.hasOption(readsOption.getOpt())) {
			READS = true;
		}
		
		SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault()
				.validationStringency(ValidationStringency.SILENT);
		
		boolean readsStartWithChr = false;
		boolean allVariantsStartWithChr = true;

		if (READS) {
			// Parse input read files and their desired min MAPQ from command line
			String[] inputREADFILEPATHS = inputREADFILESSTRING.split(",");
			String[] inputMinMAPQStrings = inputMinMAPQ.split(",");

			if (inputMinMAPQStrings.length != inputREADFILEPATHS.length) {
				throw new Exception(
						"Incorrect number of min MAPQ arguments! The number of comma separated MAPQ values must be equal to the number of provided alignment files.");
			}

			inputREADFILES = new File[inputREADFILEPATHS.length];
			minMAPQ = new double[inputMinMAPQStrings.length];
			for (int i = 0; i < inputREADFILEPATHS.length; i++) {
				minMAPQ[i] = Double.valueOf(inputMinMAPQStrings[i]);
				inputREADFILES[i] = new File(inputREADFILEPATHS[i]);
				
				if (!inputREADFILES[i].exists()) {
					throw new FileNotFoundException("File " + inputREADFILES[i].getAbsolutePath()
							+ " does not exist! All provided BAM files must be valid.");
				}
				
				SamReader samReaderCheck = samReaderFactory.open(inputREADFILES[i]);
				samReaderCheck.getFileHeader();
				SAMRecordIterator chrCheckSAMIterator = samReaderCheck.iterator();
				if(chrCheckSAMIterator.hasNext()) {
					SAMRecord firstRecord = chrCheckSAMIterator.next();
					while(firstRecord.getContig() == null) {
						if(chrCheckSAMIterator.hasNext()) {
							firstRecord = chrCheckSAMIterator.next();
						}
					}
					if(firstRecord.getContig().startsWith("chr")) {
						if(i > 0 && readsStartWithChr == false) {
							throw new Exception("All bam contigs must be uniform and start with chr or not. Mixing is not allowed.");
						}
						readsStartWithChr = true;
					} else if (readsStartWithChr == true) {
						throw new Exception("All bam contigs must be uniform and start with chr or not. Mixing is not allowed.");
					}
				} else {
					throw new Exception("Empty reads file " + inputREADFILES[i].getAbsolutePath() + ". Please provide non-empty files.");
				}
			}
			
		}

		// Ensure required files all exist
		if (!inputVCF_FILTER.exists()) {
			throw new FileNotFoundException("File " + inputVCF_FILTER.getAbsolutePath() + " does not exist!");
		}
		if (!inputVCF_ALL.exists()) {
			throw new FileNotFoundException("File " + inputVCF_ALL.getAbsolutePath() + " does not exist!");
		}
		if (!TRIO && !READS) {
			throw new FileNotFoundException(
					"You must provide at least one valid bam file containing reads (-r) or activate trio phasing (-t).");
		}

		if (OUTPUT.toFile().exists()) {
			System.out
					.println("WARNING! Output file " + OUTPUT.toString() + " already exists and will be overwritten.");
			OUTPUT.toFile().delete();
		}

		Path validationOutputPath = OUTPUT;
		Path validationOutputPathConf = OUTPUT;
		Path validationOutputPathBlockLength = OUTPUT;
		String outputPathString = OUTPUT.toString();

		if (VALIDATION) {
			int lastDot = outputPathString.lastIndexOf('.');
			validationOutputPath = Paths
					.get(outputPathString.substring(0, lastDot) + "_VALIDATION.csv");
			validationOutputPathConf = Paths.get(
					outputPathString.substring(0, lastDot) + "_VALIDATION_CONF" + outputPathString.substring(lastDot));
			validationOutputPathBlockLength = Paths.get(
					outputPathString.substring(0, lastDot) + "_VALIDATION_BLOCK_LENGTH" + outputPathString.substring(lastDot));
			if (validationOutputPath.toFile().exists()) {
				System.out.println("WARNING! Output file " + validationOutputPath.toString()
						+ " already exists and will be overwritten.");
				validationOutputPath.toFile().delete();
			}
			if (validationOutputPathConf.toFile().exists()) {
				System.out.println("WARNING! Output file " + validationOutputPathConf.toString()
						+ " already exists and will be overwritten.");
				validationOutputPathConf.toFile().delete();
			}
			if (validationOutputPathBlockLength.toFile().exists()) {
				System.out.println("WARNING! Output file " + validationOutputPathBlockLength.toString()
						+ " already exists and will be overwritten.");
				validationOutputPathBlockLength.toFile().delete();
			}
		}

		SAMFileHeader allContigsHeader = new SAMFileHeader();
		allContigsHeader.addSequence(new SAMSequenceRecord("1", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("2", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("3", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("4", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("5", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("6", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("7", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("8", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("9", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("10", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("11", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("12", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("13", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("14", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("15", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("16", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("17", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("18", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("19", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("20", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("21", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("22", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("X", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("Y", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("M", Integer.MAX_VALUE));

		IntervalList iList = new IntervalList(allContigsHeader);

		if (!PAIRED) {
			// Grab all intervals from bed file and store in interval list
			try (BufferedReader brBED = new BufferedReader(new FileReader(inputBED))) {
				String line;
				String header = null;
				String[] columns;
				int nameCol = 3;

				while ((line = brBED.readLine()) != null) {
					if (line.startsWith("#")) {
						columns = line.split("\\s");
						for (int x = 0; x < columns.length; x++) {
							if (columns[x].equals("name")) {
								nameCol = x;
							}
						}
						continue;
					}

					if (line.startsWith("track") || line.startsWith("browser") && header == null) {
						header = line;
						line = brBED.readLine();
					}
					// Read in important data from bed file and split [0] =
					// contig,
					// [1] = rangeStart, [2] = rangeEnd
					columns = line.split("\\s");

					if (columns.length < 3) {
						throw new Exception(
								"Invalid .bed file. At least three tab or white-space seperated collumns must be present.");
					}
					String name = null;
					if (columns.length > 3) {
						name = columns[nameCol];
					}
					String intervalChromosome = columns[0];
					if (intervalChromosome.startsWith("chr")) {
						intervalChromosome = intervalChromosome.substring(3, intervalChromosome.length());
					}
					iList.add(new Interval(intervalChromosome, Integer.parseInt(columns[1]),
							Integer.parseInt(columns[2]), true, name));
				}
			} catch (Exception e) {
				throw new Exception("Exception while reading bed file: \n" + e.getMessage());
			}
		}

		iList = iList.uniqued();

		// Read both VCF files
		FilteredVariantReader filteredVCFReader = new FilteredVariantReader(inputVCF_FILTER, PAIRED, PATIENT_ID, iList);
		@SuppressWarnings("resource")
		VCFFileReader allVCFReader = new VCFFileReader(inputVCF_ALL);
		VCFFileReader writeVCFReader = new VCFFileReader(inputVCF_ALL);
		CloseableIterator<VariantContext> writeVCFReadIterator = writeVCFReader.iterator();
		VariantContext curVarAllVars = writeVCFReadIterator.next();
		writeVCFReader.close();
		
		VariantContextWriter vcfWriter = null;
		String vcfOutPath = "";
		if(WRITE_VCF) {
			vcfOutPath = inputVCF_ALL.getAbsolutePath().replace(".vcf.gz", "_sp.vcf.gz");
			vcfWriter = new VariantContextWriterBuilder().setReferenceDictionary(VCFFileReader.getSequenceDictionary(inputVCF_ALL)).setOutputFile(vcfOutPath).setOutputFileType(VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF).build();
			VCFHeader header = allVCFReader.getFileHeader();
			header.addMetaDataLine(new VCFFormatHeaderLine("SPGT", 1, VCFHeaderLineType.String, "Phasing haplotype information, generated by SmartPhase, describing how the alternate alleles are phased in relation to one another"));
			header.addMetaDataLine(new VCFFormatHeaderLine("SPID", 1, VCFHeaderLineType.String, "Phasing ID information, generated by SmartPhase, where each unique ID within a given sample (but not across samples) connects records within a phasing group"));
			vcfWriter.setHeader(header);
			vcfWriter.writeHeader(header);
		}
		
		// Check if contigs start with string "chr"
		VCFFileReader chrCheckAllVCFReader = new VCFFileReader(inputVCF_ALL);
		CloseableIterator<VariantContext> chrCheckAllVCFIterator = chrCheckAllVCFReader.iterator();
		if(chrCheckAllVCFIterator.hasNext()) {
			VariantContext chrCheckVarContext = chrCheckAllVCFIterator.next();
			if(chrCheckVarContext.getContig().startsWith("chr")) {
				allVariantsStartWithChr = true;
			} else {
				allVariantsStartWithChr = false;
			}
		} else {
			chrCheckAllVCFReader.close();
			throw new Exception("All variants file must contain at least one variant!");
		}
		chrCheckAllVCFReader.close();

		if (PAIRED) {
			iList = filteredVCFReader.getiList();
			iList = iList.uniqued();
		}

		ArrayList<SamReader> samReaderSet = new ArrayList<SamReader>();
		if (READS) {
			// Create factory to read bam files and add iterators to list
			for (File inputReadsFile : inputREADFILES) {
				SamReader samReader = samReaderFactory.open(inputReadsFile);
				samReaderSet.add(samReader);
			}
		}

		ArrayList<String> confidenceValidation = new ArrayList<>();  
		ArrayList<Integer> haplotypeBlockLengths = new ArrayList<>();

		samReaderFactory = null;

		// Iterate over genomic regions
		Iterator<Interval> intervalListIterator = iList.iterator();
		ArrayList<VariantContext> variantsToPhase;

		int innocCounter = 0;

		String prevContig = "";
		Interval curInterval;
		String intervalContig;
		String intervalBamContig;
		String intervalVarContig;
		String intervalName;
		int intervalStart;
		int intervalEnd;
		String intervalIdentifier;
		boolean contigSwitch;
		
		while (intervalListIterator.hasNext()) {
			curInterval = intervalListIterator.next();
			intervalContig = curInterval.getContig();
			intervalBamContig = intervalContig;
			intervalVarContig = intervalContig;
			intervalName = "";
			intervalStart = curInterval.getStart();
			intervalEnd = curInterval.getEnd();
			intervalIdentifier = "";
			
			if(readsStartWithChr) {
				intervalBamContig = "chr"+intervalContig;
			}
			
			if(allVariantsStartWithChr) {
				intervalVarContig = "chr"+intervalContig;
			}
			
			intervalName = curInterval.getName();
			intervalIdentifier = (intervalName != null) ? intervalName + "-" + intervalContig + "-" + intervalStart + "-" + intervalEnd 
															: intervalContig + "-" + intervalStart + "-" + intervalEnd;

			if (!filteredVCFReader.contigImportantCheck(intervalContig)) {
				if(!unimportantContigs.contains(intervalContig)) {
					System.out.println("Intervals skipped because contig "+intervalContig+" not considered important. This is because the contig was either not found in the filtered variants file provided or not marked as present in the header of the VCF.");					
					unimportantContigs.add(intervalContig);
				}
				continue;
			}
			
			contigSwitch = false;
			if (!prevContig.equals(intervalContig)) { 
				contigSwitch = true; 
				prevContig = intervalContig; 
			}
			
			samIteratorList.forEach(s -> s.close());
			samIteratorList.clear();
			samIteratorList.trimToSize();
			grabLastRec.clear();
			curRecords.clear();
			curRecords.trimToSize();
			pairedEndReads.clear();
			for (SamReader sr : samReaderSet) {
				samIteratorList.add(sr.queryOverlapping(intervalBamContig, intervalStart, intervalEnd));
			}
			contigSwitch = true;

			// Grab filtered variants within current region
			ArrayList<VariantContext> regionFiltVariantList = filteredVCFReader.scan(curInterval, contigSwitch);
			if (regionFiltVariantList == null) {
				continue;
			}
			regionFiltVariantList.removeIf(v -> v.getGenotype(PATIENT_ID).isHom());

			// Ensure at least two variants in region. If not, no chance of
			// compound het. and region is removed
			if (regionFiltVariantList.size() < 2) {
				// Skip all phasing steps and move onto next interval after
				// printing into file
				for (VariantContext singleVC : regionFiltVariantList) {
					try (BufferedWriter bwOUTPUT = Files.newBufferedWriter(OUTPUT, StandardCharsets.UTF_8,
							StandardOpenOption.CREATE, StandardOpenOption.APPEND)) {

						bwOUTPUT.write(intervalIdentifier + "\t" + constructVariantString(singleVC) + "\n");

					} catch (Exception e) {
						e.printStackTrace();
						throw new IOException("Exception while writting output file: " + e.getMessage());
					}
				}
				continue;
			}

			// Grab all variants within current region
			CloseableIterator<VariantContext> regionAllVariantIterator = allVCFReader.query(intervalVarContig,
					intervalStart, intervalEnd);
			variantsToPhase = new ArrayList<VariantContext>(regionAllVariantIterator.toList());

			variantsToPhase.removeIf(v -> v.getGenotype(PATIENT_ID).isHom() || v.getGenotype(PATIENT_ID).isNoCall()
					|| ((v.getGenotype(PATIENT_ID).getAllele(0).isReference()
							|| v.getGenotype(PATIENT_ID).getAllele(0).isNoCall())
							&& (v.getGenotype(PATIENT_ID).getAllele(1).isReference()
									|| v.getGenotype(PATIENT_ID).getAllele(1).isNoCall())));
			regionAllVariantIterator.close();
			regionAllVariantIterator = null;

			System.out.println("------------------");
			System.out.println("INTERVAL CONTIG: " + intervalContig);
			System.out.println("INTERVAL START: " + intervalStart);
			System.out.println("INTERVAL END: " + intervalEnd);

			System.out.println("Filtered Variants found in interval: " + regionFiltVariantList.size());

			if (variantsToPhase.isEmpty()) {
				System.err
						.println("No variants found in VCF in this interval, but more than 2 filtered variants found.");
				continue;
			}

			neverSeenVariants.clear();
			variantsToPhase.forEach(v -> neverSeenVariants.add(v));

			ArrayList<VariantContext> trioPhasedVariants = null;

			if (TRIO) {
				trioPhasedVariants = trioPhase(variantsToPhase.iterator(), familyPed);
				System.out.println("Trio size: " + trioPhasedVariants.size());
			} else if (PHYSICAL_PHASING) {
				physicalPhasingPIDMap.clear();
				variantsToPhase.iterator().forEachRemaining(v -> fillPhysicalPhasingMap(v));
			}

			ArrayList<HaplotypeBlock> phasedVars = readPhase(variantsToPhase, curInterval, trioPhasedVariants, readsStartWithChr);
			LinkedHashSet<HaplotypeBlock> deletingDups = new LinkedHashSet<HaplotypeBlock>(phasedVars);
			phasedVars = new ArrayList<HaplotypeBlock>(deletingDups);
			for (VariantContext v : neverSeenVariants) {
				System.err.println("Never saw variant: " + v.toString());
			}

			// If trio information is available, use parents GT to resolve phase
			// where possible and then merge blocks
			if (TRIO) {
				updateTripHet(trioPhasedVariants, phasedVars);
				phasedVars = mergeBlocks(trioPhasedVariants, phasedVars);
			}

			// Write final output file
			try (final BufferedWriter bwOUTPUT = Files.newBufferedWriter(OUTPUT, StandardCharsets.UTF_8,
					StandardOpenOption.CREATE, StandardOpenOption.APPEND);
					final BufferedWriter bwVALIDATION = Files.newBufferedWriter(validationOutputPath,
							StandardCharsets.UTF_8, StandardOpenOption.CREATE, StandardOpenOption.APPEND);) {

				if (VALIDATION) {
					double allVarsToPhaseCount = (double) regionFiltVariantList.size();
					double phaseConnections = allVarsToPhaseCount - (double) phasedVars.size();
					double totalSwitchError = 0;
					double blocksExamined = 0;
					int blockLength;
					for (HaplotypeBlock phasedBlock : phasedVars) {
						blockLength = phasedBlock.getBlockEnd() - phasedBlock.getBlockStart();
						haplotypeBlockLengths.add(blockLength + 1);
						ArrayList<VariantContext> allBlockVars = phasedBlock.getAllVariants();
						if (allBlockVars.size() > 1) {
							String prevStrand = null;
							String[] prevGoldenSplit = null;
							VariantContext prevVar = null;
							double switchError = 0;
							blocksExamined++;
							allBlockVarsLoop:
							for (VariantContext var : allBlockVars) {
								VariantContext simVar = grabSimilarVarFromList(var, regionFiltVariantList);
								if(simVar != null) {
									String[] currentGoldenSplit = simVar.getGenotype(PATIENT_ID).getGenotypeString(false)
											.split("\\|");
									String currentStrand = phasedBlock.getStrand(var).toString();
									if (prevStrand != null) {
										// Check if either var belongs to notPhased group
										for (VariantContext notSeenVar : neverSeenVariants) {
											if (varsAreSimilar(notSeenVar, var) || varsAreSimilar(notSeenVar, prevVar)) {
												double confidence = phasedBlock.calculateConfidence(var, prevVar);
												if(confidence != 1) {
													phaseConnections -= 1;
													continue allBlockVarsLoop;
												}
											}
										}
										
										boolean goldenCis = false;
										
										if ((currentGoldenSplit[0].indexOf("*") != -1 && prevGoldenSplit[0].indexOf("*") != -1)
												|| (currentGoldenSplit[1].indexOf("*") != -1
												&& prevGoldenSplit[1].indexOf("*") != -1)) {
											goldenCis = true;
										}
										
										boolean predictCis = prevStrand.equals(currentStrand);
										if ((predictCis && !goldenCis) || (!predictCis && goldenCis)) {
											switchError++;
										}
									}
									prevStrand = currentStrand;
									prevGoldenSplit = currentGoldenSplit;
									prevVar = var;	
								}
							}
							totalSwitchError += switchError;
						}
					}
					double meanSwitchError = totalSwitchError/blocksExamined;
					double phasingError = meanSwitchError / phaseConnections;
					double phaseConnectionsRatio = phaseConnections / (allVarsToPhaseCount - 1.0);
					String validationString = intervalIdentifier + "," + allVarsToPhaseCount + "," + phaseConnections
							+ "," + phaseConnectionsRatio + "," + totalSwitchError + "," + meanSwitchError + "," + phasingError + "\n";
					bwVALIDATION.write(validationString);
				}

				Set<VariantContext> missingVars = new HashSet<VariantContext>();
				for (int outerCount = 0; outerCount < regionFiltVariantList.size() - 1; outerCount++) {
					VariantContext outerVariant = regionFiltVariantList.get(outerCount);
					boolean foundOuter = false;
					BitSet outerBitSet = null;
					for (int innerCount = outerCount + 1; innerCount < regionFiltVariantList.size(); innerCount++) {

						boolean foundInner = false;
						BitSet innerBitSet = null;

						boolean notPhased = true;
						boolean InnocuousFlag = false;
						boolean isTrans = false;
						boolean isCis = false;
						boolean neverSeenFlag = false;
						VariantContext innerVariant = regionFiltVariantList.get(innerCount);

						double totalConfidence = 0;
						for (HaplotypeBlock hb : phasedVars) {
							if (hb == null) {
								throw new Exception("HaplotypeBlock is null!");
							}

							VariantStrandPair varStrandPair = hb.getStrandPair(outerVariant, innerVariant);
							VariantContext trueOuterVariant = varStrandPair.getOuterVariant();
							VariantContext trueInnerVariant = varStrandPair.getInnerVariant();
							HaplotypeBlock.Strand outerStrand = varStrandPair.getOuterVariantStrand();
							HaplotypeBlock.Strand innerStrand = varStrandPair.getInnerVariantStrand();

							if (outerStrand != null) {
								foundOuter = true;
								outerBitSet = (BitSet) trueOuterVariant.getAttribute("VarFlags", null);
								if (trueOuterVariant.getAttributeAsBoolean("Innocuous", false)) {
									InnocuousFlag = true;
								}
							}

							if (innerStrand != null) {
								foundInner = true;
								innerBitSet = (BitSet) trueInnerVariant.getAttribute("VarFlags", null);
								if (trueInnerVariant.getAttributeAsBoolean("Innocuous", false)) {
									InnocuousFlag = true;
								}

							}

							if (outerStrand != null && innerStrand != null) {
								totalConfidence = hb.calculateConfidence(trueInnerVariant, trueOuterVariant);
								notPhased = false;
								isTrans = (outerStrand != innerStrand) ? true : false;
								isCis = isTrans ? false : true;

								break;
							}
						}

						// Check if inner and outer variant can be labeled as
						// innocuous based on parents GT
						if (outerBitSet != null && innerBitSet != null) {
							if ((outerBitSet.get(0) && innerBitSet.get(2)) || (outerBitSet.get(1) && innerBitSet.get(3))
									|| (innerBitSet.get(0) && outerBitSet.get(2))
									|| (innerBitSet.get(1) && outerBitSet.get(3))) {
								InnocuousFlag = true;
							}
						}


						for (VariantContext notSeenVar : neverSeenVariants) {
							if (varsAreSimilar(notSeenVar, outerVariant)
									|| varsAreSimilar(notSeenVar, innerVariant)) {
								neverSeenFlag = true;
								if (PHYSICAL_PHASING) {
									// Check if physical phasing info from GATK can
									// phase
									
									String ppInnerKey = constructPPKey(innerVariant);
									String ppOuterKey = constructPPKey(outerVariant);
									
									if (physicalPhasingPIDMap.containsKey(ppInnerKey)
											&& physicalPhasingPIDMap.containsKey(ppOuterKey) && physicalPhasingPIDMap
											.get(ppInnerKey).equals(physicalPhasingPIDMap.get(ppOuterKey))) {
										String innerVarPhase = physicalPhasingPGTMap.get(ppInnerKey);
										String outerVarPhase = physicalPhasingPGTMap.get(ppOuterKey);
										
										ppPhased++;
										
										notPhased = false;
										totalConfidence = 1;
										isTrans = (innerVarPhase.equals(outerVarPhase)) ? false : true;
										isCis = isTrans ? false : true;
									}
								}
								if (totalConfidence != 1) {
									notPhased = true;
									isTrans = false;
									isCis = false;
									totalConfidence = 0;
									break;
								}
							}
						}

						if (!foundInner) {
							missingVars.add(innerVariant);
						}

						if (InnocuousFlag) {
							innocCounter++;
						}

						// Create bitset based on booleans
						BitSet flagBits = new BitSet(5);
						flagBits.set(0, isCis);
						flagBits.set(1, isTrans);
						flagBits.set(2, notPhased);
						flagBits.set(3, InnocuousFlag);
						flagBits.set(4, neverSeenFlag);

						// Parse integer flag from bitset
						int flag = 0;
						for (int i = flagBits.nextSetBit(0); i >= 0; i = flagBits.nextSetBit(i + 1)) {
							flag += (1 << i);
						}

						// Calc total number of reads spanning both variants
						String baseOutput = intervalIdentifier + "\t" + constructVariantString(outerVariant) + "\t" + constructVariantString(innerVariant) + "\t" + flag + "\t"
								+ totalConfidence;

						if (VALIDATION) {
							String[] goldenOuterSplit = outerVariant.getGenotype(PATIENT_ID).getGenotypeString(false)
									.split("\\|");
							String[] goldenInnerSplit = innerVariant.getGenotype(PATIENT_ID).getGenotypeString(false)
									.split("\\|");

							boolean goldenPhaseCis = false;

							if ((goldenOuterSplit[0].indexOf("*") != -1 && goldenInnerSplit[0].indexOf("*") != -1)
									|| (goldenOuterSplit[1].indexOf("*") != -1
											&& goldenInnerSplit[1].indexOf("*") != -1)) {
								goldenPhaseCis = true;
							}

							if ((goldenPhaseCis && isCis) || (!goldenPhaseCis && isTrans)) {
								confidenceValidation.add(intervalIdentifier + "\t" + constructVariantString(outerVariant) + "\t" + constructVariantString(innerVariant) + "\t" + totalConfidence + "\t" + 1);
							} else if ((goldenPhaseCis && isTrans) || (!goldenPhaseCis && isCis)) {
								confidenceValidation.add(intervalIdentifier + "\t" + constructVariantString(outerVariant) + "\t" + constructVariantString(innerVariant) + "\t" + totalConfidence + "\t" + 0);
							}
						}

						if (countReads) {
							int spanningReads = countReads(curInterval, innerVariant, outerVariant);

							bwOUTPUT.write(baseOutput + "\t" + spanningReads + "\n");
						} else {
							bwOUTPUT.write(baseOutput + "\n");
						}
					}
					if (!foundOuter) {
						missingVars.add(outerVariant);
					}
				}
				
				// Write vcf file with phase
				if(WRITE_VCF) {
					while(writeVCFReadIterator.hasNext() && curVarAllVars.getStart() <= intervalEnd) {
						boolean found = false;
						for (HaplotypeBlock hb : phasedVars) {
							VariantContext phasedAllVar = hb.getSimVC(curVarAllVars);
							if(phasedAllVar != null && hb.getAllVariants().size() > 1) {
								found = true;
								if(hb.getMinConf() > vcfCutoff) {
									String phase = "";
									phase = hb.getStrand(phasedAllVar) == Strand.STRAND1 ? "0|1" : "1|0"; 
									GenotypesContext toWriteGTs = GenotypesContext.copy(curVarAllVars.getGenotypes());
									toWriteGTs.remove(toWriteGTs.get(PATIENT_ID));
									toWriteGTs.add(new GenotypeBuilder(curVarAllVars.getGenotype(PATIENT_ID)).attribute("SPGT", phase).attribute("SPID", intervalContig + "_" + hb.getBlockStart()).make());
									vcfWriter.add(new VariantContextBuilder(phasedAllVar).rmAttributes(noPrintAttributes).genotypes(toWriteGTs).make());
								}
							}	
						}
						if(!found) {
							vcfWriter.add(curVarAllVars);
						}
						curVarAllVars = writeVCFReadIterator.next();
					}
				}
				variantsToPhase = null;

				// Inform user of missing vars
				for (VariantContext missingVar : missingVars) {
					System.err.println("Could not find variant: " + constructVariantString(missingVar));
				}

			} catch (Exception e) {
				e.printStackTrace();
				throw new IOException("Exception while writting output file: " + e.getMessage());
			}

			System.out.println("------------------");
		}

		filteredVCFReader.close();
		allVCFReader.close();
		if(WRITE_VCF) {			
			vcfWriter.close();
		}

		// Time benchmarking
		long millis = System.currentTimeMillis() - startTime;
		long second = (millis / 1000) % 60;
		long minute = (millis / (1000 * 60)) % 60;
		long hour = (millis / (1000 * 60 * 60)) % 24;

		double averageCisLength = (double) globalCisLength / globalCis;
		double averageTransLength = (double) globalTransLength / globalTrans;

		// Output statistics
		String outputStatistics = "Denovo count: " + denovoCounter + "\nCis count: " + globalCis
				+ "\nAvg dist between cis: " + averageCisLength + "\nTrans count: " + globalTrans
				+ "\nAvg dist between trans: " + averageTransLength + "\nNewblock count: " + globalNewBlock
				+ "\nContradiction count: " + globalContradiction + "\nInnocuous count: " + innocCounter
				+ "\nRNAseq jump count: " + globalRNAseqCount + "\nPhysical Phasing Count: " + ppPhased + "\n"
				+ String.format("%02d:%02d:%02d:%d", hour, minute, second, millis);

		System.out.println(outputStatistics);

		if (VALIDATION) {
			try (final BufferedWriter bwVALIDATION = Files.newBufferedWriter(validationOutputPathConf,
					StandardCharsets.UTF_8, StandardOpenOption.CREATE, StandardOpenOption.APPEND)) {
				for (String conf : confidenceValidation) {
					bwVALIDATION.write(conf + "\n");
				}
			}
			try (final BufferedWriter bwVALIDATION = Files.newBufferedWriter(validationOutputPathBlockLength,
					StandardCharsets.UTF_8, StandardOpenOption.CREATE, StandardOpenOption.APPEND)) {
				bwVALIDATION.write("#Block Lengths:+\n");
				for (int blockLength : haplotypeBlockLengths) {
					bwVALIDATION.write(blockLength + "\n");
				}
			}
		}

		try (BufferedWriter bwOUTPUT = Files.newBufferedWriter(OUTPUT, StandardCharsets.UTF_8,
				StandardOpenOption.CREATE, StandardOpenOption.APPEND)) {
			PrintWriter pwOUTPUT = new PrintWriter(bwOUTPUT);
			pwOUTPUT.println(outputStatistics);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private static VariantContext grabSimilarVarFromList(VariantContext vc, ArrayList<VariantContext> regionFiltVariantList) {
		for(VariantContext posVC : regionFiltVariantList) {
			if(varsAreSimilar(vc, posVC)) {
				return posVC;
			}
		}
		return null;
	}
	
	private static boolean varsAreSimilar(VariantContext var1, VariantContext var2) {
		if(var1.getStart() == var2.getStart() && var1.getReference().equals(var2.getReference())
				&& var1.getAlternateAllele(0).equals(var2.getAlternateAllele(0))){
			return true;
		}
		return false;
	}
	
	private static String constructVariantString(VariantContext var) {
		return var.getContig() + "-" + var.getStart() + "-"
				+ var.getReference().getBaseString() + "-"
				+ var.getAlternateAllele(0).getBaseString();
	}
	
	private static String constructPPKey(VariantContext var) {
		return chrStripper(var.getContig()) + "|" + var.getStart() + "|"
				+ var.getEnd() + var.getReference().getDisplayString() + "|"
				+ var.getAlternateAllele(0).getDisplayString();
	}
	
	private static String chrStripper(String varString) {
		return varString.startsWith("chr") ? varString.substring(3) : varString;
	}

	private static int countReads(Interval interval, VariantContext vc1, VariantContext vc2) {
		int count = 0;

		if (READS) {
			SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault()
					.validationStringency(ValidationStringency.SILENT);
			SamReader samReader = null;
			int indx = 0;
			// Create factory to read bam files and add iterators to list
			for (File inputReadsFile : inputREADFILES) {
				samReader = samReaderFactory.open(inputReadsFile);
				Set<SAMRecord> srSet1 = new HashSet<SAMRecord>();
				SAMRecordIterator readIt1 = samReader.queryOverlapping(interval.getContig(), vc1.getStart(),
						vc1.getEnd());
				while (readIt1.hasNext()) {
					srSet1.add(readIt1.next());
				}
				readIt1.close();

				Set<SAMRecord> srSet2 = new HashSet<SAMRecord>();
				SAMRecordIterator readIt2 = samReader.queryOverlapping(interval.getContig(), vc2.getStart(),
						vc2.getEnd());
				while (readIt2.hasNext()) {
					srSet2.add(readIt2.next());
				}
				readIt2.close();

				srSet1.retainAll(srSet2);

				double minQ = minMAPQ[indx];
				indx++;

				Iterator<SAMRecord> readIt = srSet1.iterator();
				while (readIt.hasNext()) {
					SAMRecord curRecord = readIt.next();
					// Only use reads that are mapped and have desired quality
					if (curRecord.getReadUnmappedFlag() || curRecord.getMappingQuality() < minQ
							|| (curRecord.getReadPairedFlag() && !curRecord.getProperPairFlag())
							|| curRecord.getDuplicateReadFlag() || curRecord.isSecondaryAlignment()
							|| curRecord.getStart() > vc1.getStart() || curRecord.getEnd() < vc2.getEnd()) {
						continue;
					}
					count++;
				}
			}
		}
		return count;
	}

	public static String grabPatID() {
		return PATIENT_ID;
	}

	private static ArrayList<HaplotypeBlock> readPhase(ArrayList<VariantContext> variantsToPhase, Interval curInterval,
			ArrayList<VariantContext> trioVars, boolean readsStartWithChr) throws Exception {

		// First var is guaranteed earliest position, last var end not
		// guaranteed last, thus use interval end.
		int intervalStart = variantsToPhase.get(0).getStart();
		int intervalEnd = curInterval.getEnd();
		String intervalContig = curInterval.getContig();
		
		if(readsStartWithChr) {
			intervalContig = "chr" + intervalContig;
		}

		System.out.println("Variants found in interval: " + variantsToPhase.size());

		int readsExamined = 0;
		// Cycle through all read files
		for (int indx = 0; indx < samIteratorList.size(); indx++) {
			SAMRecordIterator curIterator = samIteratorList.get(indx);
			double minQ = minMAPQ[indx];

			// Update curRec to last record looked at by current iterator
			curRec = grabLastRec.getOrDefault(curIterator, null);

			// Add records until start of read exceeds end of interval
			while ((curRec == null || curRec.getStart() < intervalEnd) && curIterator.hasNext()) {

				curRec = curIterator.next();

				// Only use reads that are mapped and have desired quality
				if (curRec.getReadUnmappedFlag() || curRec.getMappingQuality() < minQ
						|| (curRec.getReadPairedFlag() && !curRec.getProperPairFlag()) || curRec.getDuplicateReadFlag()
						|| curRec.isSecondaryAlignment()) {
					continue;
				}
				
				if (curRec.getEnd() >= intervalStart) {
					readsExamined++;
				}

				if (curRec.getEnd() >= intervalStart) {
					if (curRec.getReadPairedFlag()) {
						pairedEndReads.put(changeNameToOtherPair(curRec.getPairedReadName()), curRec);
					}
					curRecords.add(curRec);
				}

				if (!curRec.getContig().equals(intervalContig)) {
					break;
				}
			}
			// Before switching to the next iterator, save the last record
			// looked at by current iterator
			grabLastRec.put(curIterator, curRec);
		}

		// Because variants are start sorted, reads with end less than
		// intervalStart will never be relevant again
		curRecords.removeIf(r -> r.getEnd() < intervalStart);

		// TODO: All records are saved twice now, requiring double the
		// memory..... all
		// because of one record....
		// Create temp records list to be trimmed at end.
		ArrayList<SAMRecord> trimmedRecords = new ArrayList<SAMRecord>(curRecords);

		if (curRec != null && curRec.getStart() > intervalEnd && trimmedRecords.size() > 0) {
			trimmedRecords.remove(trimmedRecords.size() - 1);
		}

		System.out.println("Reads examined in interval: " + readsExamined);
		System.out.println("Reads passing QC in interval: " + trimmedRecords.size());

		if (trimmedRecords.size() > 0) {
			return phasePIR(variantsToPhase, trimmedRecords, curInterval, trioVars);
		}

		// No reads found, so create new HB for each read so at least trio
		// phasing can be done
		ArrayList<HaplotypeBlock> intervalBlocks = new ArrayList<HaplotypeBlock>();
		for (VariantContext vc : variantsToPhase) {

			globalNewBlock++;
			// Cannot phase. Open new haplotypeBlock
			HaplotypeBlock hapBlock = new HaplotypeBlock(PATIENT_ID);
			VariantContext newVar = new VariantContextBuilder(vc)
					.genotypes(new GenotypeBuilder(vc.getGenotype(PATIENT_ID)).attribute("ReadConfidence", 1.0).make())
					.attribute("Preceding", null).make();
			hapBlock.addVariant(newVar, HaplotypeBlock.Strand.STRAND1);
			intervalBlocks.add(hapBlock);
		}
		return intervalBlocks;
	}
	
	private static String changeNameToOtherPair(String name) throws Exception {
		if(name.charAt(name.length()-3) == '1') {
			name = name.substring(0, name.length()-3) + "2" + name.substring(name.length()-2);
		} else if(name.charAt(name.length()-3) == '2') {
			name = name.substring(0, name.length()-3) + "1" + name.substring(name.length()-2);
		} else {
			throw new Exception("Expected paired name to be either 1/2 or 2/2");
		}
		return name;
	}

	private static ArrayList<HaplotypeBlock> phasePIR(ArrayList<VariantContext> variantsToPhase,
			ArrayList<SAMRecord> trimmedRecords, Interval curInterval, ArrayList<VariantContext> trioVars)
			throws Exception {

		Set<VariantContext> key;

		phaseCounter.clear();
		skipIntronCounter.clear();

		exonStartVars = new HashMap<VariantContext, VariantContext>();
		HashSet<String> examinedPairedReadNames = new HashSet<String>();

		for (SAMRecord r : trimmedRecords) {
			seenInRead.clear();
			NOT_SeenInRead.clear();
			varToStartHash.clear();
			varToEndHash.clear();
			// More than one paired end read found
			if (pairedEndReads.containsKey(r.getPairedReadName()) && pairedEndReads.get(r.getPairedReadName()) != null) {
				SAMRecord pairedRecord = pairedEndReads.get(r.getPairedReadName());
				
				//ensure paired end reads are only looked at once
				if(examinedPairedReadNames.contains(r.getReadName())) {
					continue;
				}
				
				examinedPairedReadNames.add(r.getReadName());
				
				// Grab a random variant from read to create the
				// links for merging down the road
				VariantContext ranVar1;
				VariantContext ranVar2;
				ranVar2 = countEvidence(pairedRecord, variantsToPhase, true);
				ranVar1 = countEvidence(r, variantsToPhase, true);
				if(ranVar1 != null && ranVar2 != null) {
					if(ranVar1.getStart() != ranVar2.getStart()) {
						pairedReadsVarMaps.put(ranVar1, ranVar2);
						pairedReadsVarMaps.put(ranVar2, ranVar1);
					}
				}
			} else {
				countEvidence(r, variantsToPhase, false);
			}
		}
		
		//removeNotFoundVarsFromEvidences(variantsToPhase);
		
		trimmedRecords = null;

		// No reads found spanning two vars. Creating new HB for each variant so
		// they can be trio phased
		if (phaseCounter.size() == 0) {
			// No reads found, so create new HB for each read so at least trio
			// phasing can be done
			ArrayList<HaplotypeBlock> intervalBlocks = new ArrayList<HaplotypeBlock>();
			for (VariantContext vc : variantsToPhase) {
				globalNewBlock++;
				// Cannot phase. Open new haplotypeBlock
				HaplotypeBlock hapBlock = new HaplotypeBlock(PATIENT_ID);
				VariantContext newVar = new VariantContextBuilder(vc)
						.genotypes(
								new GenotypeBuilder(vc.getGenotype(PATIENT_ID)).attribute("ReadConfidence", 1.0).make())
						.attribute("Preceding", null).make();
				hapBlock.addVariant(newVar, HaplotypeBlock.Strand.STRAND1);
				intervalBlocks.add(hapBlock);
			}
			return intervalBlocks;
		}

		double cisCounter;
		double transCounter;
		double observedCounter;
		double confidence;

		ArrayList<HaplotypeBlock> intervalBlocks = new ArrayList<HaplotypeBlock>();

		// Initialize first haplotype block and set first variant to strand1.
		// All subsequent variants phased with respect
		HaplotypeBlock hapBlock = new HaplotypeBlock(PATIENT_ID);
		VariantContext origVar = variantsToPhase.get(0);
		hapBlock.addVariant(new VariantContextBuilder(origVar)
				.genotypes(new GenotypeBuilder(origVar.getGenotype(PATIENT_ID)).attribute("ReadConfidence", 1.0).make())
				.attribute("Preceding", null).make(), HaplotypeBlock.Strand.STRAND1);
		origVar = null;

		VariantContext firstVar;
		VariantContext secondVar;

		// Examine triovars to ensure no contradictions
		VariantContext curTrioVar = null;
		VariantContext blockTrioVar1 = null;
		VariantContext blockTrioVar2 = null;

		Iterator<VariantContext> trioVarsIterator = null;
		if (TRIO) {
			trioVarsIterator = trioVars.iterator();
			if (trioVarsIterator.hasNext()) {
				curTrioVar = trioVarsIterator.next();
			}
		}

		// Phase each variant with respect to its immediate neighbor.
		for (int i = 0; i < variantsToPhase.size() - 1; i++) {
			firstVar = variantsToPhase.get(i);
			secondVar = variantsToPhase.get(i + 1);

			// Keep track if there is at least one trioVar in block
			// (blockTrioVar1) and keep an eye out
			// for any others (blockTrioVar2) to ensure no trio phasing is
			// disregarded when creating blocks
			boolean checkContradiction = false;
			if (TRIO) {
				if (blockTrioVar1 == null) {
					while (trioVarsIterator.hasNext() && curTrioVar.getStart() < firstVar.getStart()) {
						curTrioVar = trioVarsIterator.next();
					}
					if (curTrioVar != null && curTrioVar.getStart() == firstVar.getStart()
							&& curTrioVar.getGenotype(PATIENT_ID).isPhased()) {
						blockTrioVar1 = curTrioVar;
					}
				}
				if (blockTrioVar1 != null) {
					while (trioVarsIterator.hasNext() && curTrioVar.getStart() < secondVar.getStart()) {
						curTrioVar = trioVarsIterator.next();
					}
					if(firstVar.getStart() == secondVar.getStart() && trioVarsIterator.hasNext()) {
						curTrioVar = trioVarsIterator.next();
					}
					if (curTrioVar != null && curTrioVar.getStart() == secondVar.getStart()
							&& curTrioVar.getGenotype(PATIENT_ID).isPhased()) {
						blockTrioVar2 = curTrioVar;
						checkContradiction = true;
					}
				}
			}

			HaplotypeBlock.Strand firstVarStrand = hapBlock.getStrandSimVC(firstVar);
			if (firstVarStrand == null) {
				throw new Exception("firstVarStrand is null!");
			}

			key = new HashSet<VariantContext>();
			key.add(secondVar);
			key.add(firstVar);

			cisCounter = phaseCounter.getOrDefault(new PhaseCountTriple<Set<VariantContext>, Phase>(key, Phase.CIS),
					0.0);
			transCounter = phaseCounter.getOrDefault(new PhaseCountTriple<Set<VariantContext>, Phase>(key, Phase.TRANS),
					0.0);
			observedCounter = phaseCounter.getOrDefault(
					new PhaseCountTriple<Set<VariantContext>, Phase>(key, Phase.TOTAL_OBSERVED), Double.MAX_VALUE);

			confidence = Math.abs((transCounter - Math.min(2 * cisCounter, observedCounter)) / (observedCounter + 1));

			// Check if trio info contradicts how cis/trans counters would add
			// var to curBlock
			if (checkContradiction) {
				String[] prevTrioSplit = blockTrioVar1.getGenotype(PATIENT_ID).getGenotypeString(false).split("\\|");
				String[] curTrioSplit = blockTrioVar2.getGenotype(PATIENT_ID).getGenotypeString(false).split("\\|");

				double trioConf = (double) blockTrioVar1.getGenotype(PATIENT_ID).getAnyAttribute("TrioConfidence")
						* (double) blockTrioVar2.getGenotype(PATIENT_ID).getAnyAttribute("TrioConfidence");

				// CIS
				if (((prevTrioSplit[0].indexOf("*") != -1 || prevTrioSplit[0].indexOf(".") != -1) && (curTrioSplit[0].indexOf("*") != -1 || curTrioSplit[0].indexOf(".") != -1))
						|| ((prevTrioSplit[1].indexOf("*") != -1 || prevTrioSplit[1].indexOf(".") != -1) && (curTrioSplit[1].indexOf("*") != -1 || curTrioSplit[1].indexOf(".") != -1))) {
					// Contradiction
					if (transCounter >= cisCounter
							&& hapBlock.getStrandSimVC(firstVar) == hapBlock.getStrandSimVC(blockTrioVar1)) {
						// Trio has better confidence
						if (trioConf > confidence) {
							globalContradiction++;
							VariantContext newVar = new VariantContextBuilder(secondVar)
									.genotypes(new GenotypeBuilder(secondVar.getGenotype(PATIENT_ID))
											.attribute("ReadConfidence", trioConf).make())
									.attribute("Preceding", firstVar).make();

							hapBlock.addVariant(newVar, firstVarStrand);
							newVar = new VariantContextBuilder(newVar)
									.attribute("mergedBlocks", hapBlock.getHighestMB()).make();
							HaplotypeBlock returnedBlock = mergePairedReads(firstVar, secondVar, hapBlock,
									intervalBlocks, newVar);
							if (returnedBlock != null) {
								hapBlock = returnedBlock;
							}
							continue;
						}
					} else if (cisCounter >= transCounter
							&& hapBlock.getStrandSimVC(firstVar) != hapBlock.getStrandSimVC(blockTrioVar1)) {
						// Trio has better confidence
						if (trioConf > confidence) {
							globalContradiction++;
							VariantContext newVar = new VariantContextBuilder(secondVar)
									.genotypes(new GenotypeBuilder(secondVar.getGenotype(PATIENT_ID))
											.attribute("ReadConfidence", trioConf).make())
									.attribute("Preceding", firstVar).make();

							hapBlock.addVariant(newVar, hapBlock.getOppStrand(firstVarStrand));
							newVar = new VariantContextBuilder(newVar)
									.attribute("mergedBlocks", hapBlock.getHighestMB()).make();
							HaplotypeBlock returnedBlock = mergePairedReads(firstVar, secondVar, hapBlock,
									intervalBlocks, newVar);
							if (returnedBlock != null) {
								hapBlock = returnedBlock;
							}
							continue;
						}
					}
				}
				// TRANS
				else {
					// Contradiction
					if (cisCounter >= transCounter
							&& hapBlock.getStrandSimVC(firstVar) == hapBlock.getStrandSimVC(blockTrioVar1)) {
						// Trio has better confidence
						if (trioConf > confidence) {
							globalContradiction++;
							VariantContext newVar = new VariantContextBuilder(secondVar)
									.genotypes(new GenotypeBuilder(secondVar.getGenotype(PATIENT_ID))
											.attribute("ReadConfidence", trioConf).make())
									.attribute("Preceding", firstVar).make();

							hapBlock.addVariant(newVar, hapBlock.getOppStrand(firstVarStrand));
							newVar = new VariantContextBuilder(newVar)
									.attribute("mergedBlocks", hapBlock.getHighestMB()).make();
							HaplotypeBlock returnedBlock = mergePairedReads(firstVar, secondVar, hapBlock,
									intervalBlocks, newVar);
							if (returnedBlock != null) {
								hapBlock = returnedBlock;
							}
							continue;
						}
					} else if (transCounter >= cisCounter
							&& hapBlock.getStrandSimVC(firstVar) != hapBlock.getStrandSimVC(blockTrioVar1)) {
						// Trio has better confidence
						if (trioConf > confidence) {
							globalContradiction++;
							VariantContext newVar = new VariantContextBuilder(secondVar)
									.genotypes(new GenotypeBuilder(secondVar.getGenotype(PATIENT_ID))
											.attribute("ReadConfidence", trioConf).make())
									.attribute("Preceding", firstVar).make();

							hapBlock.addVariant(newVar, firstVarStrand);
							newVar = new VariantContextBuilder(newVar)
									.attribute("mergedBlocks", hapBlock.getHighestMB()).make();
							HaplotypeBlock returnedBlock = mergePairedReads(firstVar, secondVar, hapBlock,
									intervalBlocks, newVar);
							if (returnedBlock != null) {
								hapBlock = returnedBlock;
							}
							continue;
						}
					}
				}
			}

			VariantContext newVar = new VariantContextBuilder(secondVar)
					.genotypes(new GenotypeBuilder(secondVar.getGenotype(PATIENT_ID))
							.attribute("ReadConfidence", confidence).make())
					.attribute("Preceding", firstVar).make();

			if (cisCounter > transCounter) {
				globalCisLength += (secondVar.getEnd() - firstVar.getStart());
				globalCis++;
				hapBlock.addVariant(newVar, firstVarStrand);
			} else if (transCounter > cisCounter) {
				globalTransLength += (secondVar.getEnd() - firstVar.getStart());
				globalTrans++;
				hapBlock.addVariant(newVar, hapBlock.getOppStrand(firstVarStrand));
			} else {
				blockTrioVar1 = null;

				globalNewBlock++;
				// Cannot phase. Open new haplotypeBlock
				intervalBlocks.add(hapBlock);
				hapBlock = new HaplotypeBlock(PATIENT_ID);
				newVar = new VariantContextBuilder(secondVar).genotypes(
						new GenotypeBuilder(secondVar.getGenotype(PATIENT_ID)).attribute("ReadConfidence", 1.0).make())
						.attribute("Preceding", null).make();
				hapBlock.addVariant(newVar, HaplotypeBlock.Strand.STRAND1);
			}
			newVar = new VariantContextBuilder(newVar).attribute("mergedBlocks", hapBlock.getHighestMB()).make();

			// SecondVar is start of another exon. Check if can merge with other
			// block
			if (exonStartVars.containsKey(secondVar)) {
				VariantContext endOfFirstExon = exonStartVars.get(secondVar);
				key = new HashSet<VariantContext>();
				key.add(secondVar);
				key.add(endOfFirstExon);
				cisCounter = skipIntronCounter
						.getOrDefault(new PhaseCountTriple<Set<VariantContext>, Phase>(key, Phase.CIS), 0.0);
				transCounter = skipIntronCounter
						.getOrDefault(new PhaseCountTriple<Set<VariantContext>, Phase>(key, Phase.TRANS), 0.0);

				observedCounter = phaseCounter.getOrDefault(
						new PhaseCountTriple<Set<VariantContext>, Phase>(key, Phase.TOTAL_OBSERVED), Double.MAX_VALUE);

				confidence = Math
						.abs((transCounter - Math.min(2 * cisCounter, observedCounter)) / (observedCounter + 1));

				if (cisCounter != transCounter) {
					for (HaplotypeBlock hb : intervalBlocks) {
						HaplotypeBlock.Strand ex1Strand = hb.getStrandSimVC(endOfFirstExon);
						if (ex1Strand != null) {
							// Replace variant with one with a reference to link
							// for confidence calculations later
							VariantContext linkingVariant = new VariantContextBuilder(newVar)
									.attribute("linkedPreceding", hb.getSimVC(endOfFirstExon))
									.attribute("linkedConfidence", confidence).make();
							newVar = hapBlock.getSimVC(newVar);
							hapBlock.replaceVariant(linkingVariant, newVar, hapBlock.getStrand(newVar));
							globalRNAseqCount++;

							int newMergeBlock = hb.getHighestMB();
							
							newMergeBlock++;

							if (cisCounter > transCounter) {
								hb.addVariantsMerge(hapBlock.getStrandVariants(hapBlock.getStrandSimVC(secondVar)),
										ex1Strand, newMergeBlock);
								hb.addVariantsMerge(
										hapBlock.getStrandVariants(hb.getOppStrand(hapBlock.getStrandSimVC(secondVar))),
										hb.getOppStrand(ex1Strand), newMergeBlock);
								hb.updateHMC();
								hapBlock = hb;
							} else if (transCounter > cisCounter) {
								hb.addVariantsMerge(hapBlock.getStrandVariants(hapBlock.getStrandSimVC(secondVar)),
										hb.getOppStrand(ex1Strand), newMergeBlock);
								hb.addVariantsMerge(
										hapBlock.getStrandVariants(hb.getOppStrand(hapBlock.getStrandSimVC(secondVar))),
										ex1Strand, newMergeBlock);
								hb.updateHMC();
								hapBlock = hb;
							}
						}
					}
				}

			}

			HaplotypeBlock returnedBlock = mergePairedReads(firstVar, secondVar, hapBlock, intervalBlocks, newVar);
			if (returnedBlock != null) {
				hapBlock = returnedBlock;
			}

		}
		intervalBlocks.add(hapBlock);

		// TODO: Remove inefficiencies created by calculating/storing values
		// twice for overlapping intervals.
		return intervalBlocks;
	}

	// Merge blocks if variant is part of read pair
	private static HaplotypeBlock mergePairedReads(VariantContext firstVar, VariantContext secondVar,
			HaplotypeBlock hapBlock, ArrayList<HaplotypeBlock> intervalBlocks, VariantContext oldVar) throws Exception {
		ListIterator<HaplotypeBlock> intervalBlocksIterator = null;
		if (intervalBlocks.size() > 0) {
			intervalBlocksIterator = intervalBlocks.listIterator(intervalBlocks.size() - 1);
		} else {
			return null;
		}
		HaplotypeBlock hb = null;
		if (pairedReadsVarMaps.containsKey(secondVar)) {
			VariantContext connectionVar = pairedReadsVarMaps.get(secondVar);
			while (intervalBlocksIterator.hasNext()) {
				hb = intervalBlocksIterator.next();
				HaplotypeBlock.Strand conVarStrand = hb.getStrandSimVC(connectionVar);
				if (conVarStrand != null && hapBlock.getStrandSimVC(connectionVar) == null) {
					HashSet<VariantContext> key = new HashSet<VariantContext>();
					key.add(secondVar);
					key.add(connectionVar);
					double cisCounter = phaseCounter
							.getOrDefault(new PhaseCountTriple<Set<VariantContext>, Phase>(key, Phase.CIS), 0.0);
					double transCounter = phaseCounter
							.getOrDefault(new PhaseCountTriple<Set<VariantContext>, Phase>(key, Phase.TRANS), 0.0);
					double observedCounter = phaseCounter.getOrDefault(
							new PhaseCountTriple<Set<VariantContext>, Phase>(key, Phase.TOTAL_OBSERVED),
							Double.MAX_VALUE);
				
					
					double confidence = Math
							.abs((transCounter - Math.min(2 * cisCounter, observedCounter)) / (observedCounter + 1));

					// Replace variant with one with a reference to link for
					// confidence calculations later
					VariantContext linkingVariant = new VariantContextBuilder(oldVar)
							.attribute("linkedPreceding", hb.getSimVC(connectionVar))
							.attribute("linkedConfidence", confidence).make();
					int newMergeBlock = hb.getHighestMB();
					newMergeBlock++;

					oldVar = hapBlock.getSimVC(oldVar);

					if (cisCounter > transCounter) {
						intervalBlocksIterator.remove();
						hapBlock.replaceVariant(linkingVariant, oldVar, hapBlock.getStrand(oldVar));
						hb.addVariantsMerge(hapBlock.getStrandVariants(hapBlock.getStrandSimVC(secondVar)),
								conVarStrand, newMergeBlock);
						hb.addVariantsMerge(
								hapBlock.getStrandVariants(hb.getOppStrand(hapBlock.getStrandSimVC(secondVar))),
								hb.getOppStrand(conVarStrand), newMergeBlock);
						hb.updateHMC();
						return hb;
					} else if (transCounter > cisCounter) {
						intervalBlocksIterator.remove();
						hapBlock.replaceVariant(linkingVariant, oldVar, hapBlock.getStrand(oldVar));
						hb.addVariantsMerge(hapBlock.getStrandVariants(hapBlock.getStrandSimVC(secondVar)),
								hb.getOppStrand(conVarStrand), newMergeBlock);
						hb.addVariantsMerge(
								hapBlock.getStrandVariants(hb.getOppStrand(hapBlock.getStrandSimVC(secondVar))),
								conVarStrand, newMergeBlock);
						hb.updateHMC();
						return hb;
					}
				}
			}
		}
		return null;
	}

	private static VariantContext countEvidence(SAMRecord r, ArrayList<VariantContext> variantsToPhase, boolean pairedEndRead)
			throws Exception {
		// Calculate end of first exon and start of second to find vars at
		// edges for RNAseq data
		int exon1End = 0;
		int exon2Start = 0;
		VariantContext varExon1 = null;
		VariantContext varExon2 = null;
		Cigar curCigar = r.getCigar();
		if (curCigar.containsOperator(CigarOperator.N)) {
			for (CigarElement ce : curCigar.getCigarElements()) {
				if (ce.getOperator().equals(CigarOperator.N)) {
					exon2Start = exon1End + ce.getLength();
					break;
				}
				if (ce.getOperator().consumesReferenceBases()) {
					exon1End += ce.getLength();
				}
			}
		}
		ArrayList<VariantContext> trimPosVarsInRead = new ArrayList<VariantContext>(variantsToPhase);

		trimPosVarsInRead.removeIf(v -> v.getStart() < r.getAlignmentStart() || v.getEnd() > r.getAlignmentEnd());

		/*
		 * Read covers less than 2 variants and is thus not of interest. OLD. Because of
		 * paired end reads, a single variant can still technically be of use if
		 * (trimPosVarsInRead.size() < 2) { continue; }
		 */
		if (trimPosVarsInRead.size() < 1) {
			return null;
		}

		boolean varEx1Seen = false;
		VariantContext exonVar = null;
		for (VariantContext v : trimPosVarsInRead) {

			Genotype patGT = v.getGenotype(PATIENT_ID);

			// Ensure file is normalized
			if (v.getNAlleles() > 2) {
				throw new Exception("Only normalized vcf files are allowed.");
			}

			// Check if deletion is on this read and if variant is deletion
			// variant
			boolean del = (r.getReadPositionAtReferencePosition(v.getStart() + v.getReference().length() - 1,
					false) == 0) ? true : false;
			boolean delVar = (v.isSimpleDeletion()) ? true : false;

			// Determine if insert
			boolean insertVar = (v.isSimpleInsertion()) ? true : false;
			boolean insert = (r.getReadPositionAtReferencePosition(v.getStart() + v.getAlternateAllele(0).length() - 1,
					false) != r.getReadPositionAtReferencePosition(v.getStart(), false) + 1) ? true : false;

			Allele allele = patGT.getAllele(1);
			if (!allele.basesMatch(v.getAlternateAllele(0))) {
				if (patGT.getAllele(0).basesMatch(v.getAlternateAllele(0))) {
					allele = patGT.getAllele(0);
				} else {
					throw new Exception(
							"Neither patGT allele 1 nor allele 2 match the first alternative allele of the variant. Either the vcf is not normalized or this variant shouldn't be here.");
				}
			}
			// Calculate correct position in read to be compared to
			// alternative alleles in variant
			int subStrStart = r.getReadPositionAtReferencePosition(v.getStart(), false) - 1;
			int subStrEnd = r.getReadPositionAtReferencePosition(v.getStart() + allele.length() - 1, false);

			// Disregard variants who start in the middle of deletions
			// as
			// these aren't called correctly.
			if (subStrStart == -1 || subStrEnd == 0) {
				continue;
			}

			if (insert && insertVar) {
				subStrEnd = r.getReadPositionAtReferencePosition(v.getStart() + 1, false) - 1;
				if (subStrEnd == -1) {
					continue;
				}
			}
			
			varToStartHash.put(v, subStrStart);
			varToEndHash.put(v, subStrEnd);

			if (subStrStart > subStrEnd) {
				throw new Exception("Relative start of variant in read is later than end. Something went wrong.");
			}

			// Determine variants closest to intron
			if (exon1End != 0 && exon2Start != 0) {
				if (subStrStart < exon1End) {
					if (varExon1 == null) {
						varExon1 = v;
						if ((!delVar || del) && (!insertVar || insert)) {
							if (allele.basesMatch(Arrays.copyOfRange(r.getReadBases(), subStrStart, subStrEnd))) {
								neverSeenVariants.remove(v);
								varEx1Seen = true;
							}
						}
						exonVar = varExon1;
						varToStartHash.put(exonVar, subStrStart);
						varToEndHash.put(exonVar, subStrEnd);
					} else if (varExon1.getStart() < subStrStart) {
						varEx1Seen = false;
						if ((!delVar || del) && (!insertVar || insert)) {
							if (allele.basesMatch(Arrays.copyOfRange(r.getReadBases(), subStrStart, subStrEnd))) {
								neverSeenVariants.remove(v);
								varEx1Seen = true;
							}
						}
						exonVar = varExon1;
						varToStartHash.put(exonVar, subStrStart);
						varToEndHash.put(exonVar, subStrEnd);
					}
				} else if (v.getStart() > exon2Start + r.getAlignmentStart() && varExon1 != null) {
					if (varExon2 == null) {
						if ((!delVar || del) && (!insertVar || insert)) {
							if (allele.basesMatch(Arrays.copyOfRange(r.getReadBases(), subStrStart, subStrEnd))) {
								neverSeenVariants.remove(v);
								if (varEx1Seen) {
									skipIntronCounter = updatePhaseCounter(skipIntronCounter, exonVar, v, r,
											varToStartHash, varToEndHash, Phase.CIS);
								} else {
									skipIntronCounter = updatePhaseCounter(skipIntronCounter, exonVar, v, r,
											varToStartHash, varToEndHash, Phase.TRANS);
								}
							} else if (varEx1Seen) {
								skipIntronCounter = updatePhaseCounter(skipIntronCounter, exonVar, v, r, varToStartHash,
										varToEndHash, Phase.TRANS);
							}
						}
						varExon2 = v;
						exonStartVars.put(v, varExon1);
					} else if (varExon2.getStart() > v.getStart()) {
						System.err.println("START DECREASED");
					}
				}
			}

			SAMRecord dummyRead = new SAMRecord(null);
			dummyRead.setBaseQualityString("*");

			// Increase observed count
			//phaseCounter = updatePhaseCounter(phaseCounter, seenInRead, v, dummyRead, subStrStart, subStrEnd, Phase.TOTAL_OBSERVED);
			//phaseCounter = updatePhaseCounter(phaseCounter, NOT_SeenInRead, v, dummyRead, subStrStart, subStrEnd, Phase.TOTAL_OBSERVED);

			// Check alternative allele co-occurence on read
			if ((!delVar || del) && (!insertVar || insert)) {

				// Variant is found in read
				if (allele.basesMatch(Arrays.copyOfRange(r.getReadBases(), subStrStart, subStrEnd))) {
					neverSeenVariants.remove(v);

					seenInRead.add(v);
					// Increase CIS counter for all also found with this read
					//phaseCounter = updatePhaseCounter(phaseCounter, seenInRead, v, r, subStrStart, subStrEnd, Phase.CIS);

					// Increase TRANS counter for all examined and NOT
					// found on this read
					//phaseCounter = updatePhaseCounter(phaseCounter, NOT_SeenInRead, v, r, subStrStart, subStrEnd, Phase.TRANS);
				} else {

					NOT_SeenInRead.add(v);

					// Increase TRANS counter for all examined and
					// found on this read
					//phaseCounter = updatePhaseCounter(phaseCounter, seenInRead, v, r, subStrStart, subStrEnd, Phase.TRANS);
				}
			} else {
				NOT_SeenInRead.add(v);
			}
		}
		
		HashSet<VariantContext> alreadyCounted = new HashSet<VariantContext>();
		for(VariantContext seenRead : seenInRead) {
			alreadyCounted.add(seenRead);
			for(VariantContext seenRead2 : seenInRead) {
				if(!alreadyCounted.contains(seenRead2)) {
					phaseCounter = updatePhaseCounter(phaseCounter, seenRead, seenRead2, r, varToStartHash, varToEndHash, Phase.CIS);					
				}
			}
			for(VariantContext notSeenRead : NOT_SeenInRead) {
				phaseCounter = updatePhaseCounter(phaseCounter, seenRead, notSeenRead, r, varToStartHash, varToEndHash, Phase.TRANS);
			}
		}
		
		alreadyCounted.clear();
		for(VariantContext notSeenRead1 : NOT_SeenInRead) {
			alreadyCounted.add(notSeenRead1);
			for(VariantContext notSeenRead2 : NOT_SeenInRead) {
				if(!alreadyCounted.contains(notSeenRead2)) {
					phaseCounter = updatePhaseCounter(phaseCounter, notSeenRead1, notSeenRead2, r, varToStartHash, varToEndHash, Phase.TOTAL_OBSERVED);
				}
			}
		}
		
	
		if(pairedEndRead) {
			for(int index = 0; index < trimPosVarsInRead.size(); index++) {
				VariantContext possibleLinkerVar = trimPosVarsInRead.get(index);
				if(!neverSeenVariants.contains(possibleLinkerVar)) {
					return possibleLinkerVar;
				}
			}			
		}
		return null;
		
	}
	
	private static HashMap<PhaseCountTriple<Set<VariantContext>, Phase>, Double> updatePhaseCounter(
			HashMap<PhaseCountTriple<Set<VariantContext>, Phase>, Double> phaseCounter, VariantContext v1,
			VariantContext v2, SAMRecord r, HashMap<VariantContext, Integer> varToStartHash, HashMap<VariantContext, Integer> varToEndHash, Phase phase) {
		int subStrStart1 = varToStartHash.get(v1);
		int subStrStart2 = varToStartHash.get(v2);
		int subStrEnd1 = varToEndHash.get(v1);
		int subStrEnd2 = varToEndHash.get(v2);
		
		SAMRecord r1 = r;
		SAMRecord r2 = r;
		// Check if start or end is outside of record indicating paired end read that needs to be grabbed
		if(r.getStart() > v1.getStart() || v1.getStart()  > r.getEnd() || r.getStart() > v1.getStart()+(subStrEnd1-subStrStart1-1) || v1.getStart()+(subStrEnd1-subStrStart1-1) > r.getEnd()) {
			r1 = pairedEndReads.get(r.getPairedReadName());
		}
		if(r.getStart() > v2.getStart() || v2.getStart()  > r.getEnd() || r.getStart() > v2.getStart()+(subStrEnd2-subStrStart2-1) || v2.getStart()+(subStrEnd2-subStrStart2-1)> r.getEnd()) {
			r2 = pairedEndReads.get(r.getPairedReadName());
		}
		
		// One subtracted as we are working now with indexes and not substrings
		subStrEnd1 = subStrEnd1 - 1;
		subStrEnd2 = subStrEnd2 - 1;
		byte[] baseQualities1 = r1.getBaseQualities();
		byte[] baseQualities2 = r2.getBaseQualities();
		double averageQuality = 1.0;
		
		double averageQuality1 = 0.0;
		double averageQuality2 = 0.0;
		if (!r1.getBaseQualityString().equals("*")) {
			if (baseQualities1.length != r1.getReadLength()) {
				System.err.println("Base qualities array and read length not equal");
			}

			// First var
			for (int currentBase = subStrStart1; currentBase <= subStrEnd1; currentBase++) {
				double miscallProb = Math.pow(10, Double.valueOf(baseQualities1[currentBase]) / -10);
				averageQuality1 += miscallProb;
			}
			// Take average and then change from prob. error to prob. correct
			averageQuality1 = averageQuality1 / (subStrEnd1 + 1 - subStrStart1);

		}
		
		if (!r2.getBaseQualityString().equals("*")) {
			// Second var
			for (int currentBase = subStrStart2; currentBase <= subStrEnd2; currentBase++) {
				double miscallProb = Math.pow(10, Double.valueOf(baseQualities2[currentBase]) / -10);
				averageQuality2 += miscallProb;
			}
			// Take average and then change from prob. error to prob. correct
			averageQuality2 = averageQuality2 / (subStrEnd2 + 1 - subStrStart2);	
		}
		
		averageQuality = 1.0 - ((averageQuality1+averageQuality2)/2);
		final double finalAverageQuality = averageQuality;
		HashSet<VariantContext> key = new HashSet<VariantContext>();
		key.add(v1);
		key.add(v2);
		if(!phase.equals(Phase.TOTAL_OBSERVED)) {
			phaseCounter.compute(new PhaseCountTriple<Set<VariantContext>, Phase>(key, phase),
					(k, val) -> (val == null) ? finalAverageQuality : val + finalAverageQuality);			
		}
		phaseCounter.compute(new PhaseCountTriple<Set<VariantContext>, Phase>(key, Phase.TOTAL_OBSERVED),
				(k, val) -> (val == null) ? 1.0 : val + 1.0);
		return phaseCounter;
	}

	private static ArrayList<VariantContext> trioPhase(Iterator<VariantContext> inVariantsIterator, PedFile familyPed)
			throws IOException {

		// Read in .ped file and retrieve maternal/paternal IDs and see if parent is afflicted implying innoc should be turned off
		String motherID = familyPed.get(PATIENT_ID).getMaternalId();
		String fatherID = familyPed.get(PATIENT_ID).getPaternalId();
		boolean afflictedParent = true;
		if(familyPed.get(motherID) != null && familyPed.get(fatherID) != null) {
			afflictedParent = ((familyPed.get(motherID).getPhenotype().intValue() != 1) || (familyPed.get(fatherID).getPhenotype().intValue() != 1));			
		}
		

		if (motherID.equals("0") || fatherID.equals("0")) {
			System.err.println(
					"Mother and/or father not found in PED file. Either missing or PED file corrupt. Skipping trio phasing.");
			return new ArrayList<VariantContext>();
		}

		ArrayList<VariantContext> outVariants = new ArrayList<VariantContext>();
		// Iterate over all variants in region and phase the patient's genotype
		// whenever possible using trio information
		while (inVariantsIterator.hasNext()) {
			VariantContext var = inVariantsIterator.next();
			Genotype patientGT = var.getGenotype(PATIENT_ID);

			// Analyse PID and PGT
			if (PHYSICAL_PHASING) {
				fillPhysicalPhasingMap(var);
			}

			// RESULTS
			Allele motherAllele = null;
			Allele fatherAllele = null;
			double confidence = 0;

			// Read in genotypes and alleles of family
			Allele patientAllele1 = patientGT.getAllele(0);
			Allele patientAllele2 = patientGT.getAllele(1);

			Genotype motherGT = var.getGenotype(motherID);
			List<Allele> motherAlleles = motherGT.getAlleles();

			Genotype fatherGT = var.getGenotype(fatherID);
			List<Allele> fatherAlleles = fatherGT.getAlleles();

			// Father and mother must each have at least one allele. If not =>
			// denovo
			if (fatherGT.isCalled() && motherGT.isCalled()) {
				boolean father = false;
				boolean mother = false;
				if (fatherGT.countAllele(patientAllele1) < 1 && fatherGT.countAllele(patientAllele2) < 1
						&& !fatherGT.getGenotypeString().equals("*/*")) {
					father = true;
				}

				if (motherGT.countAllele(patientAllele1) < 1 && motherGT.countAllele(patientAllele2) < 1
						&& !motherGT.getGenotypeString().equals("*/*")) {
					mother = true;
				}

				if (father || mother) {
					denovoCounter++;
				}
			}

			if (motherGT.isNoCall() || fatherGT.isNoCall()) {
				continue;
			}
			
			boolean innoc = false;
			if (!afflictedParent && ((motherGT.isHomVar() || fatherGT.isHomVar())
					|| (motherGT.sameGenotype(fatherGT) && patientGT.sameGenotype(motherGT) && patientGT.isHet()))) {
				innoc = true;
			}

			boolean motherHomVar, fatherHomVar, motherContainsVar, fatherContainsVar;

			motherHomVar = (motherGT.isHomVar()) ? true : false;
			fatherHomVar = (fatherGT.isHomVar()) ? true : false;
			motherContainsVar = (motherGT.isHet() || motherGT.isHomVar()) ? true : false;
			fatherContainsVar = (fatherGT.isHet() || fatherGT.isHomVar()) ? true : false;

			// Either mother or father contain allele seen in child, but not
			// both
			if ((motherAlleles.contains(patientAllele1) && !fatherAlleles.contains(patientAllele1)
					&& fatherAlleles.contains(patientAllele2))
					|| (fatherAlleles.contains(patientAllele2) && !motherAlleles.contains(patientAllele2)
							&& motherAlleles.contains(patientAllele1))) {
				motherAllele = patientAllele1;
				fatherAllele = patientAllele2;
				confidence = 1.0;
			} else if ((motherAlleles.contains(patientAllele2) && !fatherAlleles.contains(patientAllele2)
					&& fatherAlleles.contains(patientAllele1))
					|| (fatherAlleles.contains(patientAllele1) && !motherAlleles.contains(patientAllele1)
							&& motherAlleles.contains(patientAllele2))) {
				motherAllele = patientAllele2;
				fatherAllele = patientAllele1;
				confidence = 1.0;
			}
			// Denovo mutation and other allele is found in both parents. Parent
			// that is homozygote has greater chance.
			else if (Collections.frequency(motherAlleles, patientAllele1) == 2
					&& Collections.frequency(fatherAlleles, patientAllele1) == 1) {
				motherAllele = patientAllele1;
				fatherAllele = patientAllele2;
				confidence = 0.66;
			} else if (Collections.frequency(motherAlleles, patientAllele2) == 2
					&& Collections.frequency(fatherAlleles, patientAllele2) == 1) {
				motherAllele = patientAllele2;
				fatherAllele = patientAllele1;
				confidence = 0.66;
			} else if (Collections.frequency(fatherAlleles, patientAllele1) == 2
					&& Collections.frequency(motherAlleles, patientAllele1) == 1) {
				fatherAllele = patientAllele1;
				motherAllele = patientAllele2;
				confidence = 0.66;
			} else if (Collections.frequency(fatherAlleles, patientAllele2) == 2
					&& Collections.frequency(motherAlleles, patientAllele2) == 1) {
				fatherAllele = patientAllele2;
				motherAllele = patientAllele1;
				confidence = 0.66;
			}

			BitSet varFlagBits = new BitSet(4);
			varFlagBits.set(0, motherHomVar);
			varFlagBits.set(1, fatherHomVar);
			varFlagBits.set(2, motherContainsVar);
			varFlagBits.set(3, fatherContainsVar);

			// Check if already phased
			if (!REJECT_PHASE) {
				if (patientGT.isPhased()) {
					VariantContext vc = null;
					if (innoc) {
						vc = new VariantContextBuilder(var)
								.genotypes(new GenotypeBuilder(patientGT).attribute("TrioConfidence", 1.0).make())
								.attribute("VarFlags", varFlagBits).attribute("Innocuous", true).make();
					} else {
						vc = new VariantContextBuilder(var)
								.genotypes(new GenotypeBuilder(patientGT).attribute("TrioConfidence", 1.0).make())
								.attribute("VarFlags", varFlagBits).make();

					}
					outVariants.add(vc);
					vc = null;
					continue;
				}
			}

			if (motherGT.sameGenotype(fatherGT) && patientGT.sameGenotype(motherGT) && patientGT.isHet() && !afflictedParent) {
				// Trip-het can never be phased
				outVariants.add(
						new VariantContextBuilder(var).genotypes(new GenotypeBuilder(patientGT).phased(false).make())
								.attribute("Innocuous", true).make());
				continue;
			}

			VariantContext vc;
			// Was phased
			if (motherAllele != null) {
				ArrayList<Allele> alleles = new ArrayList<Allele>();
				alleles.add(motherAllele);
				alleles.add(fatherAllele);

				Genotype phasedGT = new GenotypeBuilder(patientGT).phased(true).attribute("TrioConfidence", confidence)
						.alleles(alleles).make();

				if (innoc) {
					vc = new VariantContextBuilder(var).genotypes(phasedGT).attribute("Innocuous", true)
							.attribute("VarFlags", varFlagBits).make();
				} else {
					vc = new VariantContextBuilder(var).genotypes(phasedGT).attribute("VarFlags", varFlagBits).make();
				}

				outVariants.add(vc);
				vc = null;
			} else if (innoc) {
				vc = new VariantContextBuilder(var).attribute("Innocuous", true).attribute("VarFlags", varFlagBits)
						.make();
				outVariants.add(vc);
				vc = null;
			}

		}

		// Returns only phased variants
		return outVariants;
	}
	
	public static void fillPhysicalPhasingMap(VariantContext var) {
		Genotype patientGT = var.getGenotype(PATIENT_ID);
		if (patientGT.hasAnyAttribute("PGT") && patientGT.hasAnyAttribute("PID")) {
			String PID = (String) patientGT.getAnyAttribute("PID");
			String PGT = (String) patientGT.getAnyAttribute("PGT");
			String key = constructPPKey(var);
	
			physicalPhasingPIDMap.put(key, PID);
			physicalPhasingPGTMap.put(key, PGT);
		}
	}

	private static ArrayList<HaplotypeBlock> mergeBlocks(ArrayList<VariantContext> trioPhasedVars,
			ArrayList<HaplotypeBlock> currentBlocks) throws Exception {
		HaplotypeBlock mergeBlock = null;
		Iterator<HaplotypeBlock> hapBlockIt = currentBlocks.iterator();
		HaplotypeBlock curBlock = null;
		VariantContext prevTrioVar = null;
		int mergeBlockCntr = 2;
		int posMergeBlockCntr1;
		int posMergeBlockCntr2;

		HapBlockLoop: while (hapBlockIt.hasNext()) {
			curBlock = hapBlockIt.next();

			for (VariantContext trioVar : trioPhasedVars) {
				if (!trioVar.getGenotype(PATIENT_ID).isPhased() && trioVar.getAttributeAsBoolean("Innocuous", false)) {
					continue;
				}

				if (trioVar.getStart() > curBlock.getBlockEnd()) {
					continue HapBlockLoop;
				}

				if (curBlock.setPhased(trioVar)) {
					// Initialize mergeblock to first block containing trio var
					if (mergeBlock == null) {
						mergeBlock = curBlock;
						mergeBlockCntr = mergeBlock.getHighestMB();
						mergeBlockCntr++;
						prevTrioVar = trioVar;
						hapBlockIt.remove();
						continue HapBlockLoop;
					}
					hapBlockIt.remove();

					// [0] is always mother. [1] is always father
					String[] prevTrioSplit = prevTrioVar.getGenotype(PATIENT_ID).getGenotypeString(false).split("\\|");
					String[] curTrioSplit = trioVar.getGenotype(PATIENT_ID).getGenotypeString(false).split("\\|");

					HaplotypeBlock.Strand prevStrandMerge = mergeBlock.getStrandSimVC(prevTrioVar);
					HaplotypeBlock.Strand prevOppStrandMerge = mergeBlock.getOppStrand(prevStrandMerge);

					HaplotypeBlock.Strand strandCur = curBlock.getStrandSimVC(trioVar);
					HaplotypeBlock.Strand oppStrandCur = curBlock.getOppStrand(strandCur);

					if (prevStrandMerge == null || prevOppStrandMerge == null) {
						throw new Exception("STRAND IS NULL");
					}

					// CIS
					if (((prevTrioSplit[0].indexOf("*") != -1 || prevTrioSplit[0].indexOf(".") != -1) && (curTrioSplit[0].indexOf("*") != -1 || curTrioSplit[0].indexOf(".") != -1))
							|| ((prevTrioSplit[1].indexOf("*") != -1 || prevTrioSplit[1].indexOf(".") != -1) && (curTrioSplit[1].indexOf("*") != -1 || curTrioSplit[1].indexOf(".") != -1))) {
						posMergeBlockCntr1 = mergeBlock.addVariantsMerge(curBlock.getStrandVariants(strandCur),
								prevStrandMerge, mergeBlockCntr);
						posMergeBlockCntr2 = mergeBlock.addVariantsMerge(curBlock.getStrandVariants(oppStrandCur),
								prevOppStrandMerge, mergeBlockCntr);
						mergeBlockCntr = (posMergeBlockCntr1 > posMergeBlockCntr2) ? posMergeBlockCntr1
								: posMergeBlockCntr2;
					}
					// TRANS
					else {
						posMergeBlockCntr1 = mergeBlock.addVariantsMerge(curBlock.getStrandVariants(strandCur),
								prevOppStrandMerge, mergeBlockCntr);
						posMergeBlockCntr2 = mergeBlock.addVariantsMerge(curBlock.getStrandVariants(oppStrandCur),
								prevStrandMerge, mergeBlockCntr);
						mergeBlockCntr = (posMergeBlockCntr1 > posMergeBlockCntr2) ? posMergeBlockCntr1
								: posMergeBlockCntr2;
					}

					mergeBlockCntr++;
					prevTrioVar = trioVar;

					continue HapBlockLoop;
				}

			}
		}
		if (mergeBlock != null) {
			currentBlocks.add(mergeBlock);
		}
		return currentBlocks;
	}

	private static void updateTripHet(ArrayList<VariantContext> trioPhasedVars, ArrayList<HaplotypeBlock> currentBlocks)
			throws Exception {

		Iterator<HaplotypeBlock> hapBlockIt = currentBlocks.iterator();
		HaplotypeBlock curBlock = null;

		// Need to look at all vars in all blocks because paired-read-joined
		// blocks now can skip single blocks and all blocks aren't perfectly
		// ordered anymore
		HapBlockIterator: while (hapBlockIt.hasNext()) {
			curBlock = hapBlockIt.next();
			for (VariantContext tVar : trioPhasedVars) {
				// Slightly more efficient as those variants that have no chance
				// of being in block aren't examined
				if (tVar.getStart() > curBlock.getBlockEnd()) {
					continue HapBlockIterator;
				}
				// Check if current trio var lands in current block. If yes,
				// check
				// if triple het should be updated
				if (curBlock.getStrandSimVC(tVar) != null) {
					if (!tVar.getGenotype(PATIENT_ID).isPhased() && (tVar.getAttributeAsBoolean("Innocuous", false))) {
						curBlock.setTripHet(tVar);
					} else {
						curBlock.setPhased(tVar);
					}
				}
			}
		}
	}
}
