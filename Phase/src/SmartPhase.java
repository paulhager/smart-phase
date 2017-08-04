import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
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
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import picard.pedigree.PedFile;

public class SmartPhase {

	public enum Phase {
		CIS, TRANS, OBSERVED
	}

	// static SAMRecordIterator samIterator;
	static ArrayList<SAMRecord> curRecords = new ArrayList<SAMRecord>();
	static ArrayList<SAMRecordIterator> samIteratorList = new ArrayList<SAMRecordIterator>();
	static HashMap<SAMRecordIterator, SAMRecord> grabLastRec = new HashMap<SAMRecordIterator, SAMRecord>();
	static SAMRecord curRec = null;
	static String prevContig = "";
	static double[] minMAPQ;

	static String PATIENT_ID;

	static boolean TRIO;
	static boolean READS;

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

	public static void main(String[] args) throws Exception {
		// Parse Options
		Options options = new Options();
		options.addRequiredOption("g", "gene-regions", true,
				"Path to file containing genomic regions to be analyzed (.bed)");
		options.addRequiredOption("f", "filtered-variants", true,
				"Path to file containing patient variants filtered for significance");
		options.addRequiredOption("a", "all-variants", true, "Path to file containing all patient variants (.vcf)");
		options.addRequiredOption("p", "patient", true, "ID of patient through vcf and ped files.");
		options.addRequiredOption("o", "output", true, "Path to desired output file.");
		options.addOption("r", "reads1", true,
				"Comma seperated list of paths to files containing aligned patient reads.");
		options.addOption("m", "mapq", true,
				"Comma seperated list of mapping quality cutoff values to use when examining reads. Each value corresponds to the min MAPQ for an input BAM file.");
		options.addOption("t", false,
				"Specify if trio information is available AND contained in original-variants file provided.");
		options.addOption("d", "ped", true, "Path to file containing vcf IDs of trio.");
		options.addOption("c", false, "Is this the cohort set?");

		CommandLineParser parser = new DefaultParser();
		CommandLine cmd = parser.parse(options, args);

		long startTime = System.currentTimeMillis();

		// File inputVCF_Mot;
		// File inputVCF_Fat;

		final File inputBED = new File(cmd.getOptionValue("g"));
		final File inputVCF_FILTER = new File(cmd.getOptionValue("f"));
		final File inputVCF_ALL = new File(cmd.getOptionValue("a"));
		final String inputREADFILESSTRING = cmd.getOptionValue("r");
		final String inputMinMAPQ = cmd.getOptionValue("m");
		final File OUTPUT = new File(cmd.getOptionValue("o"));
		PATIENT_ID = cmd.getOptionValue("p");
		final boolean COHORT;
		final File inputPEDIGREE;
		PedFile familyPed = null;
		if (cmd.hasOption("t")) {
			TRIO = true;
			inputPEDIGREE = new File(cmd.getOptionValue("d"));
			familyPed = PedFile.fromFile(inputPEDIGREE, true);
		} else {
			TRIO = false;
			inputPEDIGREE = null;
		}

		if (cmd.hasOption("c")) {
			COHORT = true;
		} else {
			COHORT = false;
		}

		if (cmd.hasOption("r")) {
			READS = true;
		}
		File[] inputREADFILES = null;

		if (READS) {
			// Parse input read files and their desired min MAPQ from command
			// line
			// string
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
			}

			if (!inputREADFILES[0].exists()) {
				throw new FileNotFoundException("File " + inputREADFILES[0].getAbsolutePath()
						+ " does not exist! You must provide at least one valid bam file containing reads or activate trio phasing (-t).");
			}
		}

		// Ensure required files all exist
		if (!inputBED.exists()) {
			throw new FileNotFoundException("File " + inputBED.getAbsolutePath() + " does not exist!");
		}
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

		if (OUTPUT.exists()) {
			System.out.println(
					"WARNING! Output file " + OUTPUT.getAbsolutePath() + " already exists and will be overwritten.");
			OUTPUT.delete();
		}

		SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault()
				.validationStringency(ValidationStringency.SILENT);

		SAMFileHeader allContigsHeader = new SAMFileHeader();
		allContigsHeader.addSequence(new SAMSequenceRecord("chr1", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("chr2", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("chr3", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("chr4", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("chr5", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("chr6", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("chr7", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("chr8", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("chr9", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("chr10", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("chr11", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("chr12", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("chr13", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("chr14", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("chr15", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("chr16", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("chr17", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("chr18", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("chr19", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("chr20", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("chr21", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("chr22", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("chrX", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("chrY", Integer.MAX_VALUE));
		allContigsHeader.addSequence(new SAMSequenceRecord("chrM", Integer.MAX_VALUE));

		IntervalList iList = new IntervalList(allContigsHeader);

		if (!COHORT) {
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
					iList.add(new Interval(columns[0], Integer.parseInt(columns[1]), Integer.parseInt(columns[2]), true,
							name));
				}
			} catch (Exception e) {
				throw new Exception("Exception while reading bed file: \n" + e.getMessage());
			}
		}

		iList = iList.uniqued();

		// Read both VCF files
		FilteredVariantReader filteredVCFReader = new FilteredVariantReader(inputVCF_FILTER, COHORT, PATIENT_ID, iList);
		@SuppressWarnings("resource")
		VCFFileReader allVCFReader = new VCFFileReader(inputVCF_ALL);

		if (COHORT) {
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

		samReaderFactory = null;

		// Iterate over genomic regions
		Iterator<Interval> intervalListIterator = iList.iterator();
		ArrayList<VariantContext> variantsToPhase;
		
		int innocCounter = 0;

		String prevContig = "";
		while (intervalListIterator.hasNext()) {
			Interval curInterval = intervalListIterator.next();
			String intervalContig = curInterval.getContig();
			String intervalName = "";
			int intervalStart = curInterval.getStart();
			int intervalEnd = curInterval.getEnd();
			if (curInterval.getName() != null) {
				intervalName = curInterval.getName();
			}

			if (!filteredVCFReader.contigImportantCheck(intervalContig)) {
				continue;
			}

			boolean contigSwitch = false;
			// Contig switched.
			if (!prevContig.equals(intervalContig)) {
				for (SAMRecordIterator srIt : samIteratorList) {
					srIt.close();
				}
				samIteratorList = new ArrayList<SAMRecordIterator>();
				grabLastRec = new HashMap<SAMRecordIterator, SAMRecord>();
				curRecords = new ArrayList<SAMRecord>();
				for (SamReader sr : samReaderSet) {
					samIteratorList.add(sr.queryOverlapping(intervalContig, intervalStart, 250000000));
				}
				contigSwitch = true;
				prevContig = intervalContig;
			}

			// Grab filtered variants within current region
			ArrayList<VariantContext> regionFiltVariantList = filteredVCFReader.scan(curInterval, contigSwitch);
			if(regionFiltVariantList == null){
				continue;
			}
			regionFiltVariantList.removeIf(v -> !v.getGenotype(PATIENT_ID).isHet());

			// Ensure at least two variants in region. If not, no chance of
			// compound het. and region is removed
			if (regionFiltVariantList.size() < 2) {
				// Skip all phasing steps and move onto next interval after
				// printing into file
				for (VariantContext singleVC : regionFiltVariantList) {
					try (BufferedWriter bwOUTPUT = new BufferedWriter(new FileWriter(OUTPUT, true))) {

						bwOUTPUT.write("INTERVAL\t" + intervalContig + "\t" + intervalStart + "\t" + intervalEnd + "\t"
								+ intervalName + "\n");

						bwOUTPUT.write(singleVC.getContig() + "-" + singleVC.getStart() + "-"
								+ singleVC.getReference().getBaseString() + "-"
								+ singleVC.getAlternateAllele(0).getBaseString() + "\n");

					} catch (Exception e) {
						e.printStackTrace();
						throw new IOException("Exception while writting output file: " + e.getMessage());
					}
				}
				continue;
			}

			// Grab all variants within current region
			CloseableIterator<VariantContext> regionAllVariantIterator = allVCFReader.query(intervalContig,
					intervalStart, intervalEnd);
			variantsToPhase = new ArrayList<VariantContext>(regionAllVariantIterator.toList());
			variantsToPhase.removeIf(v -> !v.getGenotype(PATIENT_ID).isHet());
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

			ArrayList<VariantContext> trioPhasedVariants = null;

			if (TRIO) {
				trioPhasedVariants = trioPhase(variantsToPhase.iterator(), familyPed);
			}
			ArrayList<HaplotypeBlock> phasedVars = readPhase(variantsToPhase, curInterval, trioPhasedVariants);
			variantsToPhase = null;

			// If trio information is available, use parents GT to resolve phase
			// where possible and then merge blocks
			if (TRIO) {
				updateTripHet(trioPhasedVariants, phasedVars);
				phasedVars = mergeBlocks(trioPhasedVariants, phasedVars);
			}

			// Write final output file
			// Format:
			//
			// INTERVAL\t$CONTIG\t$START\t$END
			// $VAR1START-$VAR1END\t$VAR2START-$VAR2END\t$PHASEINT\t$CONFIDENCE
			//
			// $PHASEINT:
			// 0 = NOT Compound Heterozygous
			// 1 = IS Compound Heterozygous
			// -1 = No Information
			try (BufferedWriter bwOUTPUT = new BufferedWriter(new FileWriter(OUTPUT, true))) {

				bwOUTPUT.write("INTERVAL\t" + intervalContig + "\t" + intervalStart + "\t" + intervalEnd + "\t"
						+ intervalName + "\n");

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
						VariantContext innerVariant = regionFiltVariantList.get(innerCount);

						double totalConfidence = -1;
						for (HaplotypeBlock hb : phasedVars) {
							if (hb == null) {
								throw new Exception("HaplotypeBlock is null!");
							}

							// TODO: Sift through arraylists twice... very
							// inneficient
							VariantContext trueOuterVariant = hb.getSimVC(outerVariant);
							VariantContext trueInnerVariant = hb.getSimVC(innerVariant);
							HaplotypeBlock.Strand outerStrand = hb.getStrandSimVC(outerVariant);
							HaplotypeBlock.Strand innerStrand = hb.getStrandSimVC(innerVariant);

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

						if (notPhased) {
							totalConfidence = 0;
						}

						if (!foundInner) {
							missingVars.add(innerVariant);
						}
						
						if(InnocuousFlag){
							innocCounter++;
						}

						// Create bitset based on booleans
						BitSet flagBits = new BitSet(3);
						flagBits.set(0, isTrans);
						flagBits.set(1, notPhased);
						flagBits.set(2, InnocuousFlag);

						// Parse integer flag from bitset
						int flag = 0;
						for (int i = flagBits.nextSetBit(0); i >= 0; i = flagBits.nextSetBit(i + 1)) {
							flag += (1 << i);
						}

						bwOUTPUT.write(outerVariant.getContig() + "-" + outerVariant.getStart() + "-"
								+ outerVariant.getReference().getBaseString() + "-"
								+ outerVariant.getAlternateAllele(0).getBaseString() + "\t" + innerVariant.getContig()
								+ "-" + innerVariant.getStart() + "-" + innerVariant.getReference().getBaseString()
								+ "-" + innerVariant.getAlternateAllele(0).getBaseString() + "\t" + flag + "\t"
								+ totalConfidence + "\n");
					}
					if (!foundOuter) {
						missingVars.add(outerVariant);
					}
				}

				// Inform user of missing vars
				for (VariantContext missingVar : missingVars) {
					System.err.println("Could not find variant: " + missingVar.getContig() + "-" + missingVar.getStart()
							+ "-" + missingVar.getReference().getBaseString() + "-"
							+ missingVar.getAlternateAllele(0).getBaseString());
				}

			} catch (Exception e) {
				e.printStackTrace();
				throw new IOException("Exception while writting output file: " + e.getMessage());
			}

			System.out.println("------------------");
		}

		filteredVCFReader.close();
		allVCFReader.close();

		// Time benchmarking
		long millis = System.currentTimeMillis() - startTime;
		long second = (millis / 1000) % 60;
		long minute = (millis / (1000 * 60)) % 60;
		long hour = (millis / (1000 * 60 * 60)) % 24;
		
		double averageCisLength = (double) globalCisLength / globalCis;
		double averageTransLength = (double) globalTransLength / globalTrans;

		// Output statistics
		System.out.println("Denovo count: " + denovoCounter);
		System.out.println("Cis count: " + globalCis);
		System.out.println("Avg dist between cis: " + averageCisLength);
		System.out.println("Trans count: " + globalTrans);
		System.out.println("Avg dist between trans: " + averageTransLength);
		System.out.println("Newblock count: " + globalNewBlock);
		System.out.println("Contradiction count: " + globalContradiction);
		System.out.println("Innocuous count: "+innocCounter);
		System.out.println(String.format("%02d:%02d:%02d:%d", hour, minute, second, millis));

		try (BufferedWriter bwOUTPUT = new BufferedWriter(new FileWriter(OUTPUT, true))) {
			PrintWriter pwOUTPUT = new PrintWriter(bwOUTPUT);
			pwOUTPUT.println("Denovo count: " + denovoCounter);
			pwOUTPUT.println("Cis count: " + globalCis);
			pwOUTPUT.println("Avg dist between cis: " + averageCisLength);
			pwOUTPUT.println("Trans count: " + globalTrans);
			pwOUTPUT.println("Avg dist between trans: " + averageTransLength);
			pwOUTPUT.println("Newblock count: " + globalNewBlock);
			pwOUTPUT.println("Contradiction count: " + globalContradiction);
			pwOUTPUT.println("Innocuous count: "+innocCounter);
			pwOUTPUT.println(String.format("%02d:%02d:%02d:%d", hour, minute, second, millis));
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static String grabPatID() {
		return PATIENT_ID;
	}

	private static ArrayList<HaplotypeBlock> readPhase(ArrayList<VariantContext> variantsToPhase, Interval curInterval,
			ArrayList<VariantContext> trioVars) throws Exception {

		// First var is guaranteed earliest position, last var end not
		// guaranteed last, thus use interval end.
		int intervalStart = variantsToPhase.get(0).getStart();
		int intervalEnd = curInterval.getEnd();
		String intervalContig = curInterval.getContig();

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

				if (curRec.getEnd() >= intervalStart) {
					readsExamined++;
				}

				// Only use reads that are mapped and have desired quality
				if (curRec.getReadUnmappedFlag() || curRec.getMappingQuality() < minQ
						|| (curRec.getReadPairedFlag() && !curRec.getProperPairFlag()) || curRec.getDuplicateReadFlag()
						|| curRec.getNotPrimaryAlignmentFlag()) {
					continue;
				}

				if (curRec.getEnd() >= intervalStart) {
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

	private static ArrayList<HaplotypeBlock> phasePIR(ArrayList<VariantContext> variantsToPhase,
			ArrayList<SAMRecord> trimmedRecords, Interval curInterval, ArrayList<VariantContext> trioVars)
			throws Exception {

		ArrayList<VariantContext> trimPosVarsInRead;
		Set<VariantContext> key;

		HashMap<PhaseCountTriple<Set<VariantContext>, Phase>, Integer> phaseCounter = new HashMap<PhaseCountTriple<Set<VariantContext>, Phase>, Integer>();
		HashMap<PhaseCountTriple<Set<VariantContext>, Phase>, Integer> skipIntronCounter = new HashMap<PhaseCountTriple<Set<VariantContext>, Phase>, Integer>();
		
		HashMap<VariantContext, VariantContext> exonStartVars = new HashMap<VariantContext, VariantContext>();
		
		for (SAMRecord r : trimmedRecords) {
			
			// Calculate end of first exon and start of second to find vars at edges for RNAseq data
			int exon1End = 0;
			int exon2Start = 0;
			VariantContext varExon1 = null;
			VariantContext varExon2 = null;
			Cigar curCigar = r.getCigar();
			if(curCigar.containsOperator(CigarOperator.N)){
				for(CigarElement ce : curCigar.getCigarElements()){
					if(ce.getOperator().equals(CigarOperator.N)){
						exon2Start = exon1End + ce.getLength();
						break;
					}
					if(ce.getOperator().consumesReferenceBases()){
						exon1End += ce.getLength();
					}
				}
			}
			trimPosVarsInRead = new ArrayList<VariantContext>(variantsToPhase);

			trimPosVarsInRead.removeIf(v -> v.getStart() < r.getAlignmentStart() || v.getEnd() > r.getAlignmentEnd());

			// Read covers less than 2 variants and is thus not of interest.
			if (trimPosVarsInRead.size() < 2) {
				continue;
			}
			
			
			if(exon1End != 0 && exon2Start != 0 && r.getAlignmentStart() == 256249){
				System.out.println("Alignment Start: "+r.getAlignmentStart());
				System.out.println("Exon 1 End: "+(r.getAlignmentStart()+exon1End));
				System.out.println("Exon 2 Start: "+(r.getAlignmentStart()+exon2Start));
				System.out.println("Alignment End: "+r.getAlignmentEnd());
				for (VariantContext v : trimPosVarsInRead) {
					System.out.println(v.toStringDecodeGenotypes());
				}
				System.out.println("----");
			}
			

			ArrayList<VariantContext> seenInRead = new ArrayList<VariantContext>();
			ArrayList<VariantContext> NOT_SeenInRead = new ArrayList<VariantContext>();

			
			boolean varEx1Seen = false;
			ArrayList<VariantContext> exVarList = new ArrayList<VariantContext>();
			for (VariantContext v : trimPosVarsInRead) {

				Genotype patGT = v.getGenotype(PATIENT_ID);

				// Ensure file is normalized
				if (v.getNAlleles() > 2) {
					throw new Exception("Only normalized vcf files are allowed.");
				}

				// Check if deletion is on this read and if variant is deletion
				// variant
				boolean del = (r.getReadPositionAtReferencePosition(v.getEnd(), false) == 0) ? true : false;
				boolean delVar = (v.isSimpleDeletion()) ? true : false;

				// Determine if insert
				boolean insertVar = (v.isSimpleInsertion()) ? true : false;
				boolean insert = (r.getReadPositionAtReferencePosition(v.getStart() + 1,
						false) != r.getReadPositionAtReferencePosition(v.getStart(), false) + 1) ? true : false;

				Allele allele = patGT.getAllele(1);
				if (!allele.isNonReference()) {
					throw new Error("Only normalized vcf files are allowed.");
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
				
				// Determine variants closest to intron
				if(exon1End != 0 && exon2Start != 0){
					System.out.println(r.getAlignmentStart());
					System.out.println(r.getAlignmentEnd());
					System.out.println(subStrStart);
					System.out.println(exon1End);
					System.out.println(exon2Start);
					if(subStrStart < exon1End){
						if(varExon1 == null){
							System.out.println("Set varExon1!");
							varExon1 = v;
							if ((!delVar || del) && (!insertVar || insert)) {
								if (allele.basesMatch(Arrays.copyOfRange(r.getReadBases(), subStrStart, subStrEnd))) {
									varEx1Seen = true;
								}
							}
							exVarList.clear();
							exVarList.add(varExon1);
						} else if(varExon1.getStart() < subStrStart){
							varEx1Seen = false;
							if ((!delVar || del) && (!insertVar || insert)) {
								if (allele.basesMatch(Arrays.copyOfRange(r.getReadBases(), subStrStart, subStrEnd))) {
									varEx1Seen = true;
								}
							}
							varExon1 = v;
							exVarList.clear();
							exVarList.add(varExon1);
						}
					} else if (subStrStart > exon2Start && varExon1 != null){
						if(varExon2 == null){
							System.out.println("Set varExon2!");
							if ((!delVar || del) && (!insertVar || insert)) {
								if (allele.basesMatch(Arrays.copyOfRange(r.getReadBases(), subStrStart, subStrEnd))) {
									if(varEx1Seen){
										skipIntronCounter = updatePhaseCounter(skipIntronCounter, exVarList, v, Phase.CIS);
									} else {
										skipIntronCounter = updatePhaseCounter(skipIntronCounter, exVarList, v, Phase.TRANS);
									}
								} else if (varEx1Seen){
									skipIntronCounter = updatePhaseCounter(skipIntronCounter, exVarList, v, Phase.TRANS);
								}
							}
							varExon2 = v;
							exonStartVars.put(v, varExon1);
						} else if(varExon2.getStart() > subStrStart){
							System.err.println("START DECREASED");
						}
					}
					System.out.println("---");
				}
				

				// Check alternative allele co-occurence on read
				if ((!delVar || del) && (!insertVar || insert)) {

					// Increase observed count
					phaseCounter = updatePhaseCounter(phaseCounter, seenInRead, v, Phase.OBSERVED);
					phaseCounter = updatePhaseCounter(phaseCounter, NOT_SeenInRead, v, Phase.OBSERVED);

					// Variant is found in read
					if (allele.basesMatch(Arrays.copyOfRange(r.getReadBases(), subStrStart, subStrEnd))) {

						seenInRead.add(v);
						// Increase CIS counter for all also found with this
						// read
						phaseCounter = updatePhaseCounter(phaseCounter, seenInRead, v, Phase.CIS);

						// Increase TRANS counter for all examined and NOT
						// found on this read
						phaseCounter = updatePhaseCounter(phaseCounter, NOT_SeenInRead, v, Phase.TRANS);
					} else {

						NOT_SeenInRead.add(v);

						// Increase TRANS counter for all examined and
						// found on this read
						phaseCounter = updatePhaseCounter(phaseCounter, seenInRead, v, Phase.TRANS);
					}

				}
			}
		}
		trimPosVarsInRead = null;
		trimmedRecords = null;

		// No reads found spanning two vars. Creating new HB for each variant so
		// thex can be trio phased
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
		VariantContext firstTrioVar = null;
		VariantContext secondTrioVar = null;

		Iterator<VariantContext> trioVarsIterator = null;
		if (TRIO) {
			trioVarsIterator = trioVars.iterator();
			if (trioVarsIterator.hasNext()) {
				firstTrioVar = trioVarsIterator.next();
			}
		}

		// Phase each variant with respect to its immediate neighbor.
		for (int i = 0; i < variantsToPhase.size() - 1; i++) {
			firstVar = variantsToPhase.get(i);
			secondVar = variantsToPhase.get(i + 1);

			// Check if current two vars directly next to each other whose read
			// counts are being compared are also trio vars
			boolean checkContradiction = false;
			if (TRIO) {
				if (firstTrioVar != null) {
					while ((firstTrioVar.getStart() < firstVar.getStart()
							|| !firstTrioVar.getGenotype(PATIENT_ID).isPhased()) && trioVarsIterator.hasNext()) {
						firstTrioVar = trioVarsIterator.next();
					}
					if (firstTrioVar.getStart() == firstVar.getStart() && trioVarsIterator.hasNext()) {
						secondTrioVar = trioVarsIterator.next();
						while (trioVarsIterator.hasNext() && !secondTrioVar.getGenotype(PATIENT_ID).isPhased()) {
							secondTrioVar = trioVarsIterator.next();
						}
						if (secondTrioVar.getStart() == secondVar.getStart()
								&& secondTrioVar.getGenotype(PATIENT_ID).isPhased()
								&& firstTrioVar.getGenotype(PATIENT_ID).isPhased()) {
							checkContradiction = true;
						} else {
							firstTrioVar = secondTrioVar;
							secondTrioVar = null;
						}
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

			cisCounter = phaseCounter.getOrDefault(new PhaseCountTriple<Set<VariantContext>, Phase>(key, Phase.CIS), 0);
			transCounter = phaseCounter.getOrDefault(new PhaseCountTriple<Set<VariantContext>, Phase>(key, Phase.TRANS),
					0);
			observedCounter = phaseCounter.getOrDefault(
					new PhaseCountTriple<Set<VariantContext>, Phase>(key, Phase.OBSERVED), Integer.MAX_VALUE);

			key = new HashSet<VariantContext>();
			double confidence = Math.abs((transCounter - Math.min(2*cisCounter, observedCounter)) / (observedCounter + 1));
			
			// Check if trio info contradicts cis/trans counters
			if (checkContradiction) {
				String[] prevTrioSplit = firstTrioVar.getGenotype(PATIENT_ID).getGenotypeString(false).split("\\|");
				String[] curTrioSplit = secondTrioVar.getGenotype(PATIENT_ID).getGenotypeString(false).split("\\|");

				double trioConf = (double) firstTrioVar.getGenotype(PATIENT_ID).getAnyAttribute("TrioConfidence")
						* (double) secondTrioVar.getGenotype(PATIENT_ID).getAnyAttribute("TrioConfidence");

				// CIS
				if ((prevTrioSplit[0].indexOf("*") != -1 && curTrioSplit[0].indexOf("*") != -1)
						|| (prevTrioSplit[1].indexOf("*") != -1 && curTrioSplit[1].indexOf("*") != -1)) {
					// Contradiction
					if (transCounter > cisCounter) {
						// Trio has better confidence
						if (trioConf > confidence) {
							globalContradiction++;
							VariantContext newVar = new VariantContextBuilder(secondVar)
									.genotypes(new GenotypeBuilder(secondVar.getGenotype(PATIENT_ID))
											.attribute("ReadConfidence", trioConf).make())
									.attribute("Preceding", firstVar).make();

							hapBlock.addVariant(newVar, firstVarStrand);
							continue;
						}
					}
				}
				// TRANS
				else {
					// Contradiction
					if (cisCounter > transCounter) {
						// Trio has better confidence
						if (trioConf > confidence) {
							globalContradiction++;
							VariantContext newVar = new VariantContextBuilder(secondVar)
									.genotypes(new GenotypeBuilder(secondVar.getGenotype(PATIENT_ID))
											.attribute("ReadConfidence", trioConf).make())
									.attribute("Preceding", firstVar).make();

							hapBlock.addVariant(newVar, hapBlock.getOppStrand(firstVarStrand));
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
				globalNewBlock++;
				// Cannot phase. Open new haplotypeBlock
				intervalBlocks.add(hapBlock);
				hapBlock = new HaplotypeBlock(PATIENT_ID);
				VariantContext newVarNewBlock = new VariantContextBuilder(secondVar).genotypes(
						new GenotypeBuilder(secondVar.getGenotype(PATIENT_ID)).attribute("ReadConfidence", 1.0).make())
						.attribute("Preceding", null).make();
				hapBlock.addVariant(newVarNewBlock, HaplotypeBlock.Strand.STRAND1);
			}
			
			// SecondVar is start of another exon. Check if can merge with other block
			if(exonStartVars.containsKey(secondVar)){
				VariantContext endOfFirstExon = exonStartVars.get(secondVar);
				key.add(secondVar);
				key.add(endOfFirstExon);
				cisCounter = skipIntronCounter.getOrDefault(new PhaseCountTriple<Set<VariantContext>, Phase>(key, Phase.CIS), 0);
				transCounter = skipIntronCounter.getOrDefault(new PhaseCountTriple<Set<VariantContext>, Phase>(key, Phase.TRANS),
						0);
				System.out.println("Found second key");
				if(cisCounter>transCounter){
					for(HaplotypeBlock hb : intervalBlocks){
						HaplotypeBlock.Strand ex1Strand = hb.getStrand(endOfFirstExon);
						if(ex1Strand != null){
							int origIndex = intervalBlocks.indexOf(hb);
							hb.addVariantsMerge(hapBlock.getStrandVariants(hapBlock.getStrand(secondVar)), ex1Strand, -2);
							hb.addVariantsMerge(hapBlock.getStrandVariants(hb.getOppStrand(hapBlock.getStrand(secondVar))), hb.getOppStrand(ex1Strand), -2);
							intervalBlocks.set(origIndex, hb);
						}
					}
				} else {
					for(HaplotypeBlock hb : intervalBlocks){
						HaplotypeBlock.Strand ex1Strand = hb.getStrand(endOfFirstExon);
						if(ex1Strand != null){
							int origIndex = intervalBlocks.indexOf(hb);
							hb.addVariantsMerge(hapBlock.getStrandVariants(hapBlock.getStrand(secondVar)), hb.getOppStrand(ex1Strand), -2);
							hb.addVariantsMerge(hapBlock.getStrandVariants(hb.getOppStrand(hapBlock.getStrand(secondVar))), ex1Strand, -2);
							intervalBlocks.set(origIndex, hb);
						}
					}
				}
				
			}
		}
		intervalBlocks.add(hapBlock);

		// TODO: Remove inefficiencies created by calculating/storing values
		// twice for overlapping intervals.
		return intervalBlocks;
	}

	private static HashMap<PhaseCountTriple<Set<VariantContext>, Phase>, Integer> updatePhaseCounter(
			HashMap<PhaseCountTriple<Set<VariantContext>, Phase>, Integer> phaseCounter,
			ArrayList<VariantContext> group, VariantContext v, Phase phase) {
		for (VariantContext member : group) {
			HashSet<VariantContext> key = new HashSet<VariantContext>();
			key.add(v);
			key.add(member);
			phaseCounter.compute(new PhaseCountTriple<Set<VariantContext>, Phase>(key, phase),
					(k, val) -> (val == null) ? 1 : val + 1);
		}
		return phaseCounter;
	}

	private static ArrayList<VariantContext> trioPhase(Iterator<VariantContext> inVariantsIterator, PedFile familyPed)
			throws IOException {
		// Read in .ped file and retrieve maternal/paternal IDs
		String motherID = familyPed.get(PATIENT_ID).getMaternalId();
		String fatherID = familyPed.get(PATIENT_ID).getPaternalId();

		ArrayList<VariantContext> outVariants = new ArrayList<VariantContext>();
		// Iterate over all variants in region and phase the patient's genotype
		// whenever possible using trio information
		while (inVariantsIterator.hasNext()) {
			VariantContext var = inVariantsIterator.next();
			Genotype patientGT = var.getGenotype(PATIENT_ID);

			// Check if already phased
			if (patientGT.isPhased()) {
				VariantContext vc = new VariantContextBuilder(var)
						.genotypes(new GenotypeBuilder(patientGT).attribute("TrioConfidence", 1.0).make()).make();
				outVariants.add(vc);
				vc = null;
				continue;
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
					/*
					 * pw.println("---"); pw.println("Patient: " +
					 * patientGT.toString()); pw.println("Mother: " +
					 * motherGT.toString()); pw.println("Father: " +
					 * fatherGT.toString()); pw.print("Alleles: "); for(Allele a
					 * : var.getAlleles()){ pw.print(a.getDisplayString()+" ");
					 * } pw.println(); if(father){ pw.print("Father: "); for(int
					 * pl : fatherGT.getAD()){ pw.print(pl+" "); } pw.println();
					 * }
					 * 
					 * if(mother){ pw.print("Mother: "); for(int pl :
					 * motherGT.getAD()){ pw.print(pl+" "); } pw.println(); }
					 * 
					 * System.out.println("---"); System.out.println("Patient: "
					 * + patientGT.toString()); System.out.println("Mother: " +
					 * motherGT.toString()); System.out.println("Father: " +
					 * fatherGT.toString()); System.out.print("Alleles: ");
					 * for(Allele a : var.getAlleles()){
					 * System.out.print(a.getDisplayString()+" "); }
					 * System.out.println();
					 * 
					 * if(father){ System.out.print("Father: "); boolean
					 * simpleCall = true; for(int pl : fatherGT.getAD()){ if(pl
					 * != 0){ if(!simpleCall){ notSimpleCounter++; } simpleCall
					 * = false; } System.out.print(pl+" "); }
					 * System.out.println(); }
					 * 
					 * if(mother){ boolean simpleCall = true;
					 * System.out.print("Mother: "); for(int pl :
					 * motherGT.getAD()){ if(pl != 0){ if(!simpleCall){
					 * notSimpleCounter++; } simpleCall = false; }
					 * System.out.print(pl+" "); } System.out.println(); }
					 * 
					 * 
					 */
					denovoCounter++;
				}
			}

			if (motherGT.isNoCall() || fatherGT.isNoCall()) {
				continue;
			}

			if (motherGT.sameGenotype(fatherGT) && patientGT.sameGenotype(motherGT)) {
				outVariants.add(new VariantContextBuilder(var).attribute("Innocuous", true).make());
				continue;
			}

			boolean innoc = false;
			if ((motherGT.isHomVar() && fatherGT.isHet()) || (fatherGT.isHomVar() && motherGT.isHet())
					|| (motherGT.isHomVar() && fatherGT.isHomVar())) {
				innoc = true;
			}

			boolean motherHomVar, fatherHomVar, motherContainsVar, fatherContainsVar;

			motherHomVar = (motherGT.isHomVar()) ? true : false;
			fatherHomVar = (fatherGT.isHomVar()) ? true : false;
			motherContainsVar = (motherGT.isHet() || motherGT.isHomVar()) ? true : false;
			fatherContainsVar = (fatherGT.isHet() || fatherGT.isHomVar()) ? true : false;

			motherGT = null;
			fatherGT = null;

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

	private static ArrayList<HaplotypeBlock> mergeBlocks(ArrayList<VariantContext> trioPhasedVars,
			ArrayList<HaplotypeBlock> currentBlocks) throws Exception {
		HaplotypeBlock mergeBlock = null;
		Iterator<HaplotypeBlock> hapBlockIt = currentBlocks.iterator();
		HaplotypeBlock curBlock = null;
		VariantContext prevTrioVar = null;
		if (hapBlockIt.hasNext()) {
			curBlock = hapBlockIt.next();
		} else {
			return currentBlocks;
		}

		int mergeBlockCntr = 2;
		for (VariantContext trioVar : trioPhasedVars) {

			if (!trioVar.getGenotype(PATIENT_ID).isPhased() && trioVar.getAttributeAsBoolean("Innocuous", false)) {
				continue;
			}

			// Increment blocks as long as var is ahead of block
			while (trioVar.getStart() > curBlock.getBlockEnd() && hapBlockIt.hasNext()) {
				curBlock = hapBlockIt.next();
			}
			// Check if current trio var lands in current block. If yes, merge
			while (curBlock.setPhased(trioVar) && hapBlockIt.hasNext()) {

				// Initialize mergeblock to first block containing trio var
				if (mergeBlock == null) {
					mergeBlock = curBlock;
					prevTrioVar = trioVar;
					hapBlockIt.remove();
					curBlock = hapBlockIt.next();
					continue;
				}

				// [0] is always mother. [1] is always father
				String[] prevTrioSplit = prevTrioVar.getGenotype(PATIENT_ID).getGenotypeString(false).split("\\|");
				String[] curTrioSplit = trioVar.getGenotype(PATIENT_ID).getGenotypeString(false).split("\\|");

				HaplotypeBlock.Strand prevStrand = mergeBlock.getStrandSimVC(prevTrioVar);
				HaplotypeBlock.Strand prevOppStrand = mergeBlock.getOppStrand(prevStrand);

				if (prevStrand == null || prevOppStrand == null) {
					throw new Exception("STRAND IS NULL");
				}

				// CIS
				if ((prevTrioSplit[0].indexOf("*") != -1 && curTrioSplit[0].indexOf("*") != -1)
						|| (prevTrioSplit[1].indexOf("*") != -1 && curTrioSplit[1].indexOf("*") != -1)) {
					mergeBlock.addVariantsMerge(curBlock.getStrandVariants(prevStrand), prevStrand, mergeBlockCntr);
					mergeBlock.addVariantsMerge(curBlock.getStrandVariants(prevOppStrand), prevOppStrand,
							mergeBlockCntr);
				}
				// TRANS
				else {
					mergeBlock.addVariantsMerge(curBlock.getStrandVariants(prevStrand), prevOppStrand, mergeBlockCntr);
					mergeBlock.addVariantsMerge(curBlock.getStrandVariants(prevOppStrand), prevStrand, mergeBlockCntr);
				}

				mergeBlockCntr++;
				prevTrioVar = trioVar;
				hapBlockIt.remove();
				curBlock = hapBlockIt.next();
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

		if (hapBlockIt.hasNext()) {
			curBlock = hapBlockIt.next();
		} else {
			return;
		}

		for (VariantContext trioVar : trioPhasedVars) {
			// Increment blocks as long as var is ahead of block
			while (trioVar.getStart() > curBlock.getBlockEnd() && hapBlockIt.hasNext()) {
				curBlock = hapBlockIt.next();
			}

			// Check if current trio var lands in current block. If yes, check
			// if triple het should be updated
			if (curBlock.getStrandSimVC(trioVar) != null) {
				if (!trioVar.getGenotype(PATIENT_ID).isPhased()
						&& (trioVar.getAttributeAsBoolean("Innocuous", false))) {
					curBlock.setTripHet(trioVar);
				}
			}
		}
	}

}
