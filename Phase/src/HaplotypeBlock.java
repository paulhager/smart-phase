import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class HaplotypeBlock {

	public enum Strand {
		STRAND1, STRAND2
	}

	private ArrayList<VariantContext> strand1;
	private ArrayList<VariantContext> strand2;
	private HashMap<Strand, ArrayList<VariantContext>> strandVariants;
	private int blockStart = Integer.MAX_VALUE;
	private int blockEnd = Integer.MIN_VALUE;
	private int highestMergedBlockCounter = 1;
	private int possibleHMBC = Integer.MIN_VALUE;
	private String PATIENT_ID;

	public HaplotypeBlock(String patID) {
		strand1 = new ArrayList<VariantContext>();
		strand2 = new ArrayList<VariantContext>();
		strandVariants = new HashMap<Strand, ArrayList<VariantContext>>();
		strandVariants.put(Strand.STRAND1, strand1);
		strandVariants.put(Strand.STRAND2, strand2);
		PATIENT_ID = patID;
	}

	/**
	 * Adds variant to specified strand and update blockStart or blockEnd as
	 * needed
	 * 
	 * @param vc
	 *            - Variant to add
	 * @param s
	 *            - Strand to add to
	 */
	public void addVariant(VariantContext vc, Strand s) {
		vc = new VariantContextBuilder(vc).attribute("mergedBlocks", highestMergedBlockCounter).make();
		strandVariants.get(s).add(vc);

		if (vc.getStart() < blockStart) {
			blockStart = vc.getStart();
		}

		if (vc.getEnd() > blockEnd) {
			blockEnd = vc.getEnd();
		}
	}

	/**
	 * Adds variant to specified strand and replace existing variant
	 * 
	 * @param vc
	 *            - Variant to add
	 * @param toBeReplaced
	 *            - Variant to be replaced
	 * @param s
	 *            - Strand to add to
	 */
	public void replaceVariant(VariantContext vc, VariantContext toBeReplaced, Strand s) {
		strandVariants.get(s).set(strandVariants.get(s).indexOf(toBeReplaced), vc);
	}

	/**
	 * Add multiple variants to strand and adjusts mergeBlock
	 * 
	 * @param vcList
	 *            - List of variant contexts to add
	 * @param s
	 *            - Strand to be added to
	 * @throws Exception
	 */
	public int addVariantsMerge(ArrayList<VariantContext> vcList, Strand s, int mergeBlockCntr) throws Exception {
		if (mergeBlockCntr > highestMergedBlockCounter) {
			highestMergedBlockCounter = mergeBlockCntr;
		}

		int highestBC = highestMergedBlockCounter;

		for (VariantContext vc : vcList) {
			int oldVCMergedBlock = 1;
			if (vc.hasAttribute("mergedBlocks")) {
				oldVCMergedBlock = (int) vc.getAttribute("mergedBlocks");
			} else {
				throw new Exception("Variant context is missing mergedBlocks attribute!");
			}

			// -1 because mergedBlocks start at 1
			int counterToInsert = highestMergedBlockCounter - 1 + oldVCMergedBlock;
			if (counterToInsert > highestBC) {
				highestBC = counterToInsert;
			}

			vc = new VariantContextBuilder(vc).attribute("mergedBlocks", counterToInsert).make();
			strandVariants.get(s).add(vc);

			if (vc.getStart() < blockStart) {
				blockStart = vc.getStart();
			}

			if (vc.getEnd() > blockEnd) {
				blockEnd = vc.getEnd();
			}
		}
		possibleHMBC = highestBC;
		return highestBC;
	}

	/**
	 * Get all variants
	 * 
	 * @return allVars - Concatenated arrayList of strand1 + strand2
	 */
	public ArrayList<VariantContext> getAllVariants() {
		ArrayList<VariantContext> allVars = new ArrayList<VariantContext>(strand1);
		allVars.addAll(strand2);
		return allVars;
	}

	/**
	 * Returns strand1 or strand2 based on desired strand
	 * 
	 * @param s
	 *            - Strand enum
	 * @return strandVariants.get(s)
	 */
	public ArrayList<VariantContext> getStrandVariants(Strand s) {
		return strandVariants.get(s);
	}

	/**
	 * Returns strand on which variant lays or null if variant not in block
	 * 
	 * @param vc
	 *            - Variant context for which the strand name should be returned
	 * @return Strand
	 */
	public Strand getStrand(VariantContext vc) {
		if (strand1.contains(vc)) {
			return Strand.STRAND1;
		} else if (strand2.contains(vc)) {
			return Strand.STRAND2;
		} else {
			return null;
		}
	}

	/**
	 * To be used when searching for a variant context that has been trimmed and
	 * thus the exact object can not be found (i.e. by using getStrand)
	 * 
	 * @param vc
	 *            - Variant context for which the strand name should be returned
	 * @return Strand
	 */
	public Strand getStrandSimVC(VariantContext vc) {

		for (VariantContext posVC : strand1) {
			if (posVC.getStart() == vc.getStart() && posVC.getContig().equals(vc.getContig())
					&& posVC.getReference().equals(vc.getReference())
					&& posVC.getAlternateAllele(0).equals(vc.getAlternateAllele(0))) {
				return Strand.STRAND1;
			}
		}
		for (VariantContext posVC : strand2) {
			if (posVC.getStart() == vc.getStart() && posVC.getContig().equals(vc.getContig())
					&& posVC.getReference().equals(vc.getReference())
					&& posVC.getAlternateAllele(0).equals(vc.getAlternateAllele(0))) {
				return Strand.STRAND2;
			}
		}
		return null;
	}

	/**
	 * Grabs saved variantContext with just start position
	 * 
	 * @param vc
	 *            - VariantContext with same start position as the one in block
	 *            to be grabbed
	 * @return VariantContext
	 */
	public VariantContext getSimVC(VariantContext vc) {
		for (VariantContext posVC : strand1) {
			if (posVC.getStart() == vc.getStart() && posVC.getContig().equals(vc.getContig())
					&& posVC.getReference().equals(vc.getReference())
					&& posVC.getAlternateAllele(0).equals(vc.getAlternateAllele(0))) {
				return posVC;
			}
		}
		for (VariantContext posVC : strand2) {
			if (posVC.getStart() == vc.getStart() && posVC.getContig().equals(vc.getContig())
					&& posVC.getReference().equals(vc.getReference())
					&& posVC.getAlternateAllele(0).equals(vc.getAlternateAllele(0))) {
				return posVC;
			}
		}
		return null;
	}

	/**
	 * Returns start position of earliest variant
	 * 
	 * @return blockStart
	 */
	public int getBlockStart() {
		return blockStart;
	}

	/**
	 * Returns farthest end position of a variant
	 * 
	 * @return blockEnd
	 */
	public int getBlockEnd() {
		return blockEnd;
	}

	/**
	 * Returns opposite strand
	 * 
	 * @return Strand != inputStrand
	 */
	public Strand getOppStrand(Strand s) {
		if (s == Strand.STRAND1) {
			return Strand.STRAND2;
		} else {
			return Strand.STRAND1;
		}
	}

	/**
	 * Calculates confidence for the output phase of two variants.
	 * 
	 * @param vc1
	 * @param vc2
	 * @return totalConfidence
	 * @throws Exception
	 *             mergedBlocks attribute not found
	 */
	public double calculateConfidence(VariantContext vc1, VariantContext vc2) throws Exception {
		VariantContext nearestTrioVC1 = findNearestTPhased(vc1);
		VariantContext nearestTrioVC2 = findNearestTPhased(vc2);

		// No trio info available
		if (nearestTrioVC1 == null || nearestTrioVC2 == null) {
			return multiplyConfidence(vc1, vc2).confidence();
		}

		// Get final product confidence score using trio
		ConfidencePair<Double, Integer> cpTrio1 = multiplyConfidence(vc1, nearestTrioVC1);
		cpTrio1.setConfidence(cpTrio1.confidence()
				* (double) nearestTrioVC1.getGenotype(PATIENT_ID).getAnyAttribute("TrioConfidence"));

		ConfidencePair<Double, Integer> cpTrio2 = multiplyConfidence(vc2, nearestTrioVC2);
		cpTrio1.setConfidence(cpTrio2.confidence()
				* (double) nearestTrioVC2.getGenotype(PATIENT_ID).getAnyAttribute("TrioConfidence"));

		// Try to calculate confidence using only read conf. scores if both
		// variants are in same mergeBlock.
		if (vc1.getAttributeAsInt("mergedBlocks", -1) == vc2.getAttributeAsInt("mergedBlocks", -1)) {
			ConfidencePair<Double, Integer> cpRead = multiplyConfidence(vc1, vc2);

			// Compare amount of steps to determine which conf. score is
			// returned
			if (cpRead.steps() < (cpTrio1.steps() + cpTrio2.steps()) / 2) {
				return cpRead.confidence();
			}
		}

		return cpTrio1.confidence() * cpTrio2.confidence();
	}

	private VariantContext findNearestTPhased(VariantContext vc) throws Exception {
		VariantContext nearestTrioVC = null;
		// Find nearest trioPhase within block
		if (vc.getGenotype(PATIENT_ID).isPhased()) {
			nearestTrioVC = vc;
		} else {
			int mergeB = vc.getAttributeAsInt("mergedBlocks", -1);
			if (mergeB == -1) {
				throw new Exception("mergedBlocks NOT found! Something went wrong.");
			}
			int distance = Integer.MAX_VALUE;
			for (VariantContext curVC : this.getAllVariants()) {
				if (curVC.getGenotype(PATIENT_ID).isPhased() && curVC.getAttributeAsInt("mergedBlocks", -1) == mergeB
						&& distance > Math.abs(vc.getStart() - curVC.getStart())) {
					distance = Math.abs(vc.getStart() - curVC.getStart());
					nearestTrioVC = curVC;
				}
			}
		}

		return nearestTrioVC;
	}

	private ConfidencePair<Double, Integer> multiplyConfidence(VariantContext vc1, VariantContext vc2)
			throws Exception {
		if (vc1 == vc2) {
			return new ConfidencePair<Double, Integer>(1.0, 0);
		}

		// Variants must be in same mergeBlock, or else must be linker
		if (vc1.getAttributeAsInt("mergedBlocks", -1) != vc2.getAttributeAsInt("mergedBlocks", -1)) {
			return calcConfidenceLinkedBlockPairs(vc1, vc2);
		}

		return calcConfCleanRunLinkedPreceding(vc1, vc2, null);
	}

	private ConfidencePair<Double, Integer> calcConfidenceLinkedBlockPairs(VariantContext var1, VariantContext var2)
			throws Exception {
		int var1MB = var1.getAttributeAsInt("mergedBlocks", 0);
		int var2MB = var2.getAttributeAsInt("mergedBlocks", 0);

		int higherMB;
		int lowerMB;

		boolean foundRun1 = false;
		boolean foundRun2 = false;

		ConfidencePair<Double, Integer> confP1 = null;
		ConfidencePair<Double, Integer> confP2 = null;

		int curMergedBlock;
		int linkedPrecedingCount = 1;
		ArrayList<VariantContext> tempVarConList = new ArrayList<VariantContext>();
		ArrayList<VariantContext> tempTrioVarList = new ArrayList<VariantContext>();
		ArrayList<Integer[]> tempBlockConnectionsList = new ArrayList<Integer[]>();

		HashMap<Integer, ArrayList<VariantContext>> mergedBlockVars = new HashMap<Integer, ArrayList<VariantContext>>();
		// Run through all vars and group by mergedBlock
		int hMBC = Integer.MIN_VALUE;
		for (VariantContext var : this.getAllVariants()) {
			curMergedBlock = var.getAttributeAsInt("mergedBlocks", 0);
			if (curMergedBlock > hMBC) {
				hMBC = curMergedBlock;
			}
			// Grabs the arraylist of vars with curMergedBlock or creates
			tempVarConList = (mergedBlockVars.containsKey(curMergedBlock)) ? mergedBlockVars.get(curMergedBlock)
					: new ArrayList<VariantContext>();
			tempVarConList.add(var);
			mergedBlockVars.put(curMergedBlock, tempVarConList);
		}
		// Start walking backwards through all merged block sequences and notate
		// their connections and trio vars they contain
		for (int curMergeBC = hMBC; curMergeBC > 0; curMergeBC--) {
			tempVarConList = mergedBlockVars.get(curMergeBC);
			if (curMergeBC == var1MB) {
				foundRun1 = true;
			} else if (curMergeBC == var2MB) {
				foundRun2 = true;
			}
			// Search for linkedPairs and trioVars
			for (VariantContext blockVar : tempVarConList) {
				if (blockVar.getGenotype(PATIENT_ID).hasAnyAttribute("TrioConfidence")) {
					tempTrioVarList.add(blockVar);
				}
				if (blockVar.hasAttribute("linkedPreceding")) {
					// [0] = Upstream MB, [1] = Downstream MB
					Integer downstreamMB = this.getSimVC((VariantContext) blockVar.getAttribute("linkedPreceding"))
							.getAttributeAsInt("mergedBlocks", 0);
					Integer[] curMergeBCLPTuple = { curMergeBC, downstreamMB};
					if(!tempBlockConnectionsList.contains(curMergeBCLPTuple)){
						linkedPrecedingCount++;
						tempBlockConnectionsList.add(curMergeBCLPTuple);
					}
				}
			}
			linkedPrecedingCount--;
			// End of run
						if (linkedPrecedingCount == 0) {
							linkedPrecedingCount = 1;
							// Both vars found in run, calc distance between
							if (foundRun1 && foundRun2) {
								if (var1MB > var2MB) {
									higherMB = var1MB;
									lowerMB = var2MB;
								} else {
									higherMB = var2MB;
									lowerMB = var1MB;
								}
								// Check if easy or hard case
								boolean unbrokenRun = true;
								for (int runMBCheck = higherMB; runMBCheck > lowerMB; runMBCheck--) {
									Integer[] testArray = { runMBCheck, runMBCheck - 1 };
									if (!deepContains(tempBlockConnectionsList, testArray)) {
										unbrokenRun = false;
									}
								}
								// Normal (easy) case
								if (unbrokenRun) {
									return calcConfCleanRunLinkedPreceding(var1, var2, null);
								} else {
									// There are at least two branches of LP and the vars
									// are on different branches (difficult case)

									// Find nearest block that both vars connect to
									Set<Integer> var1LinkedBlocks = new HashSet<Integer>();
									Set<Integer> var2LinkedBlocks = new HashSet<Integer>();
									Set<Integer> intersectionLinkedBlocks = new HashSet<Integer>();
									var1LinkedBlocks.add(var1MB);
									var2LinkedBlocks.add(var2MB);
									// Create set of blocks that var1 goes through
									var1LinkedBlocks = findLinkedBlocks(var1MB, tempBlockConnectionsList, var1LinkedBlocks);
									var2LinkedBlocks = findLinkedBlocks(var2MB, tempBlockConnectionsList, var2LinkedBlocks);
									// intersectionLinkedBlocks now contains intersection
									// blocks
									intersectionLinkedBlocks = new HashSet<>(var1LinkedBlocks);
									intersectionLinkedBlocks.retainAll(var2LinkedBlocks);
									// determine closest block
									int closestCommonBlock = Integer.MAX_VALUE;
									for (int blockNumber : intersectionLinkedBlocks) {
										if (blockNumber < closestCommonBlock) {
											closestCommonBlock = blockNumber;
										}
									}
									// Grab variant from closest block
									VariantContext commonBlockVar1 = null;
									VariantContext commonBlockVar2 = null;
									if(closestCommonBlock > 4){
										System.out.println();
									}
									for (VariantContext closestCommonVars : mergedBlockVars.get(closestCommonBlock)) {
										if (closestCommonVars.hasAttribute("linkedPreceding")) {
											int connectedBlock = this
													.getSimVC((VariantContext) closestCommonVars.getAttribute("linkedPreceding"))
													.getAttributeAsInt("mergedBlocks", 0);
											if (var1LinkedBlocks.contains(connectedBlock) || var1MB == closestCommonBlock) {
												commonBlockVar1 = closestCommonVars;
											}
											if (var2LinkedBlocks.contains(connectedBlock) || var2MB == closestCommonBlock) {
												commonBlockVar2 = closestCommonVars;
											}
										}
									}
									if (commonBlockVar1 == null || commonBlockVar2 == null) {
										throw new Exception("commonBlockVar1 or commonBlockVar2 is null!");
									}

									// Calc confidence between var and its connecting var in
									// the common block
									confP1 = calcConfCleanRunLinkedPreceding(var1, commonBlockVar1, var1LinkedBlocks);
									confP2 = calcConfCleanRunLinkedPreceding(var2, commonBlockVar2, var2LinkedBlocks);

									double finalConf = confP1.confidence() * confP2.confidence();
									int finalSteps = confP1.steps() + confP2.steps();

									return new ConfidencePair<Double, Integer>(finalConf, finalSteps);
								}
							} else if (foundRun1) {
								VariantContext earliestTrioVar = null;
								int earliestMergeBlock = Integer.MAX_VALUE;
								for (VariantContext posEarliestTrio : tempTrioVarList) {
									int posEarliestMB = posEarliestTrio.getAttributeAsInt("mergedBlocks", Integer.MAX_VALUE);
									if (posEarliestMB < earliestMergeBlock) {
										earliestMergeBlock = posEarliestMB;
										earliestTrioVar = posEarliestTrio;
									}
								}
								confP1 = calcConfCleanRunLinkedPreceding(var1, earliestTrioVar, null);
							} else if (foundRun2) {
								VariantContext earliestTrioVar = null;
								int earliestMergeBlock = Integer.MAX_VALUE;
								for (VariantContext posEarliestTrio : tempTrioVarList) {
									int posEarliestMB = posEarliestTrio.getAttributeAsInt("mergedBlocks", Integer.MAX_VALUE);
									if (posEarliestMB < earliestMergeBlock) {
										earliestMergeBlock = posEarliestMB;
										earliestTrioVar = posEarliestTrio;
									}
								}
								confP2 = calcConfCleanRunLinkedPreceding(var2, earliestTrioVar, null);
							}
							tempTrioVarList = new ArrayList<VariantContext>();
							tempBlockConnectionsList = new ArrayList<Integer[]>();
							foundRun1 = false;
							foundRun2 = false;
						}
		}
		if (confP1 == null || confP2 == null) {
			throw new Exception("Either both var MBs weren't found, or necessary trio vars missing");
		}

		double finalConf = confP1.confidence() * confP2.confidence();
		int finalSteps = confP1.steps() + confP2.steps();

		return new ConfidencePair<Double, Integer>(finalConf, finalSteps);
	}
	
	public static boolean deepContains(List<Integer[]> list, Integer[] probe) {
	    for (Integer[] element : list) {
	      if (Arrays.deepEquals(element, probe)) {
	        return true;
	      }
	    }
	    return false;
	  }

	private Set<Integer> findLinkedBlocks(int searchedForBlock, ArrayList<Integer[]> tempBlockConnectionsList,
			Set<Integer> resultSet) {
		Set<Integer> tempResults = new HashSet<Integer>();
		for (Integer[] linkedPair : tempBlockConnectionsList) {
			if (linkedPair[1] == searchedForBlock) {
				int connectedBlock = linkedPair[0];
				if(!resultSet.contains(connectedBlock)){
					resultSet.add(connectedBlock);
					tempResults = findLinkedBlocks(connectedBlock, tempBlockConnectionsList, resultSet);
					resultSet.addAll(tempResults);					
				}
			} else if(linkedPair[0] == searchedForBlock){
				int connectedBlock = linkedPair[1];
				if(!resultSet.contains(connectedBlock)){
					resultSet.add(connectedBlock);
					tempResults = findLinkedBlocks(connectedBlock, tempBlockConnectionsList, resultSet);
					resultSet.addAll(tempResults);					
				}
			}
		}
		return resultSet;
	}

	private ConfidencePair<Double, Integer> calcConfCleanRunLinkedPreceding(VariantContext vc1, VariantContext vc2, Set<Integer> varLinkedBlocks)
			throws Exception {
		double product = 1.0;

		if (vc1.getStart() == vc2.getStart() && vc1.getContig().equals(vc2.getContig())
				&& vc1.getReference().equals(vc2.getReference())
				&& vc1.getAlternateAllele(0).equals(vc2.getAlternateAllele(0))) {
			return new ConfidencePair<Double, Integer>(1.0, 0);
		}

		// Variants must be in same mergeBlock, or else must be linker
		if (vc1.getAttributeAsInt("mergedBlocks", -1) != vc2.getAttributeAsInt("mergedBlocks", -1)) {
			// Handle funny jumps due to paired-reads or RNAseq
			VariantContext linkerBlockVC = null;
			VariantContext otherBlockVC = null;
			if (vc1.getAttributeAsInt("mergedBlocks", -2) > vc2.getAttributeAsInt("mergedBlocks", -1)) {
				linkerBlockVC = vc1;
				otherBlockVC = vc2;
			} else if (vc2.getAttributeAsInt("mergedBlocks", -2) > vc1.getAttributeAsInt("mergedBlocks", -1)) {
				linkerBlockVC = vc2;
				otherBlockVC = vc1;
			} else {
				throw new Exception("Expecting linker var exception: \n" + vc1.toStringDecodeGenotypes() + "\n"
						+ vc2.toStringDecodeGenotypes());
			}

			VariantContext linkerV = null;;
			if(varLinkedBlocks != null){
				linkerV = calcComplicatedLinkerVar(linkerBlockVC, varLinkedBlocks);
			} else {
				linkerV = calcLinkerVar(linkerBlockVC);
			}

			ConfidencePair<Double, Integer> confP1 = calcConfCleanRunLinkedPreceding(linkerV, linkerBlockVC, varLinkedBlocks);
			ConfidencePair<Double, Integer> confP2 = calcConfCleanRunLinkedPreceding(otherBlockVC,
					getSimVC((VariantContext) linkerV.getAttribute("linkedPreceding")), varLinkedBlocks);

			double finalConf = confP1.confidence() * confP2.confidence()
					* linkerV.getAttributeAsDouble("linkedConfidence", -1.0);
			int finalSteps = 1 + confP1.steps() + confP2.steps();

			return new ConfidencePair<Double, Integer>(finalConf, finalSteps);
		}

		int vc1Start = vc1.getStart();
		int vc2Start = vc2.getStart();
		int end;
		VariantContext curVC;

		// curVC is vc with higher start, end with lower.
		end = (vc1Start > vc2Start) ? vc2Start : vc1Start;
		curVC = (vc1Start > vc2Start) ? vc1 : vc2;

		// Go from furthest VC to earliest VC, multiplying confidence along the
		// way
		int cnt = 0;
		while (curVC.getStart() != end) {
			product = product * (double) curVC.getGenotype(PATIENT_ID).getAnyAttribute("ReadConfidence");
			cnt++;
			curVC = this.getSimVC((VariantContext) curVC.getAttribute("Preceding"));
		}
		product = product * (double) curVC.getGenotype(PATIENT_ID).getAnyAttribute("ReadConfidence");

		return new ConfidencePair<Double, Integer>(product, cnt);
	}

	private VariantContext calcComplicatedLinkerVar(VariantContext linkerBlockVC, Set<Integer> varLinkedBlocks) throws Exception {
		for (VariantContext v : this.getAllVariants()) {
			if (v.hasAttribute("linkedPreceding")
					&& v.getAttributeAsInt("mergedBlocks", -3) == linkerBlockVC.getAttributeAsInt("mergedBlocks", -4)
						&& varLinkedBlocks.contains(getSimVC(((VariantContext) v.getAttribute("linkedPreceding", null))).getAttributeAsInt("mergedBlocks", -1))) {
				return v;
			}
		}
		throw new Exception(
				"NO LINKER VAR FOUND. Searching in block of: " + linkerBlockVC.toStringDecodeGenotypes());

	}

	private VariantContext calcLinkerVar(VariantContext linkerBlockVC) throws Exception {
		//ArrayList<VariantContext> linkerVar = new ArrayList<VariantContext>();
		/*
		if(linkerBlockVC.hasAttribute("linkedPreceding")){
			return linkerBlockVC;
		}
		*/

		for (VariantContext v : this.getAllVariants()) {
			if (v.hasAttribute("linkedPreceding")
					&& v.getAttributeAsInt("mergedBlocks", -3) == linkerBlockVC.getAttributeAsInt("mergedBlocks", -4)
						&& getSimVC(((VariantContext) v.getAttribute("linkedPreceding", null))).getAttributeAsInt("mergedBlocks", -1) == (linkerBlockVC.getAttributeAsInt("mergedBlocks", -2)-1)) {
				return v;
			}
		}
		throw new Exception(
					"NO LINKER VAR FOUND. Searching in block of: " + linkerBlockVC.toStringDecodeGenotypes());
		/*
		if (linkerVar.size() == 1){
			return linkerVar.get(0);
		}
		VariantContext downStreamLV = null;
		int dSLVStart = Integer.MAX_VALUE;
		for(VariantContext posLinkerVar : linkerVar){
			if(posLinkerVar.getStart() < dSLVStart){
				dSLVStart = posLinkerVar.getStart();
				downStreamLV = posLinkerVar;
			}
		}
		return downStreamLV;
		*/
	}

	/**
	 * Set phased attribute of a variant's patient_GT context to true.
	 * 
	 * @param trioVar
	 *            - Variant whose phase is to be changed
	 * @return True - If var found False - If var not found
	 */
	public boolean setPhased(VariantContext trioVar) {
		for (VariantContext posVC : strand1) {
			if (posVC.getStart() == trioVar.getStart()) {
				if (trioVar.getAttributeAsBoolean("Innocuous", false)) {
					this.replaceVariant(new VariantContextBuilder(posVC)
							.genotypes(new GenotypeBuilder(posVC.getGenotype(PATIENT_ID)).phased(true)
									.attribute("TrioConfidence",
											trioVar.getGenotype(PATIENT_ID).getAnyAttribute("TrioConfidence"))
									.make())
							.attribute("Innocuous", true).make(), posVC, Strand.STRAND1);
				} else {
					this.replaceVariant(new VariantContextBuilder(posVC)
							.genotypes(new GenotypeBuilder(posVC.getGenotype(PATIENT_ID)).phased(true)
									.attribute("TrioConfidence",
											trioVar.getGenotype(PATIENT_ID).getAnyAttribute("TrioConfidence"))
									.make())
							.make(), posVC, Strand.STRAND1);
				}
				return true;
			}
		}
		for (VariantContext posVC : strand2) {
			if (posVC.getStart() == trioVar.getStart()) {
				if (trioVar.getAttributeAsBoolean("Innocuous", false)) {
					this.replaceVariant(new VariantContextBuilder(posVC)
							.genotypes(new GenotypeBuilder(posVC.getGenotype(PATIENT_ID)).phased(true)
									.attribute("TrioConfidence",
											trioVar.getGenotype(PATIENT_ID).getAnyAttribute("TrioConfidence"))
									.make())
							.attribute("Innocuous", true).make(), posVC, Strand.STRAND2);
				} else {
					this.replaceVariant(new VariantContextBuilder(posVC)
							.genotypes(new GenotypeBuilder(posVC.getGenotype(PATIENT_ID)).phased(true)
									.attribute("TrioConfidence",
											trioVar.getGenotype(PATIENT_ID).getAnyAttribute("TrioConfidence"))
									.make())
							.make(), posVC, Strand.STRAND2);
				}
				return true;
			}
		}
		return false;
	}

	/**
	 * Set Innocuous attribute of a variant's patient_GT context to true.
	 * 
	 * @param trioVar
	 *            - Variant whose phase is to be changed
	 * @return True - If var found False - If var not found
	 */
	public boolean setTripHet(VariantContext trioVar) {
		for (VariantContext posVC : strand1) {
			if (posVC.getStart() == trioVar.getStart()) {
				this.replaceVariant(new VariantContextBuilder(posVC).attribute("Innocuous", true).make(), posVC,
						Strand.STRAND1);
				return true;
			}
		}
		for (VariantContext posVC : strand2) {
			if (posVC.getStart() == trioVar.getStart()) {
				this.replaceVariant(new VariantContextBuilder(posVC).attribute("Innocuous", true).make(), posVC,
						Strand.STRAND2);
				return true;
			}
		}
		return false;
	}

	/**
	 * Returns highestMergedBlockCounter value.
	 * 
	 * @return highestMergedBlockCounter
	 */
	public int getHighestMB() {
		return highestMergedBlockCounter;
	}
	
	public void updateHMC(){
		if(possibleHMBC > highestMergedBlockCounter){
			highestMergedBlockCounter = possibleHMBC;
		}
	}

}
