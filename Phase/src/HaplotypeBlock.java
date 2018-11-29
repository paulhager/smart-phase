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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

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
			if (posVC.getStart() == vc.getStart() && posVC.getReference().equals(vc.getReference())
					&& posVC.getAlternateAllele(0).equals(vc.getAlternateAllele(0))) {
				return Strand.STRAND1;
			}
		}
		for (VariantContext posVC : strand2) {
			if (posVC.getStart() == vc.getStart() && posVC.getReference().equals(vc.getReference())
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
			if (posVC.getStart() == vc.getStart() && posVC.getReference().equals(vc.getReference())
					&& posVC.getAlternateAllele(0).equals(vc.getAlternateAllele(0))) {
				return posVC;
			}
		}
		for (VariantContext posVC : strand2) {
			if (posVC.getStart() == vc.getStart() && posVC.getReference().equals(vc.getReference())
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
		if(nearestTrioVC1.getGenotype(PATIENT_ID) == null || cpTrio1 == null || nearestTrioVC1 == null || nearestTrioVC1.getGenotype(PATIENT_ID).getAnyAttribute("TrioConfidence") == null){
			System.out.println("hi");
		}
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
		if (vc.getGenotype(PATIENT_ID).isPhased() && vc.getGenotype(PATIENT_ID).hasAnyAttribute("TrioConfidence")) {
			nearestTrioVC = vc;
		} else {
			int mergeB = vc.getAttributeAsInt("mergedBlocks", -1);
			if (mergeB == -1) {
				throw new Exception("mergedBlocks NOT found! Something went wrong.");
			}
			int distance = Integer.MAX_VALUE;
			for (VariantContext curVC : this.getAllVariants()) {
				if (curVC.getGenotype(PATIENT_ID).isPhased() && curVC.getGenotype(PATIENT_ID).hasAnyAttribute("TrioConfidence") && curVC.getAttributeAsInt("mergedBlocks", -1) == mergeB
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

		return calcConfCleanRunLinkedPreceding(vc1, vc2, 0);
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
					Integer[] curMergeBCLPTuple = { curMergeBC, downstreamMB };
					if (!tempBlockConnectionsList.contains(curMergeBCLPTuple)) {
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
						return calcConfCleanRunLinkedPreceding(var1, var2, 0);
					} else {
						// There are at least two branches of LP and the vars
						// are on different branches (difficult case)

						// Determine connecting path
						ArrayList<Integer> connectingPath = calculateConnectingPath(var1MB, var2MB,
								tempBlockConnectionsList, new ArrayList<Integer>());

						ConfidencePair<Double, Integer> resultConfidence = new ConfidencePair<Double, Integer>(1.0, 0);
						VariantContext varInMBVar1 = null;
						VariantContext varInMBVar2 = null;

						// Calculate confidence between each block in path
						Integer block1 = connectingPath.get(0);
						for (int step = 1; step < connectingPath.size(); step++) {
							Integer block2 = connectingPath.get(step);
							Integer upstreamBlock = 0;
							Integer downstreamBlock = 0;
							if (block2 > block1) {
								upstreamBlock = block2;
								downstreamBlock = block1;
							} else {
								upstreamBlock = block1;
								downstreamBlock = block2;
							}
							for (VariantContext posConVar : mergedBlockVars.get(upstreamBlock)) {
								if (posConVar.hasAttribute("linkedPreceding")) {
									VariantContext linkedVar = this
											.getSimVC((VariantContext) posConVar.getAttribute("linkedPreceding"));
									int linkedBlock = linkedVar.getAttributeAsInt("mergedBlocks", -1);
									if (linkedBlock == downstreamBlock) {
										// Grab first and last var for final
										// connections
										if (step == 1 && linkedBlock == block1) {
											varInMBVar1 = linkedVar;
										} else if (step == 1) {
											varInMBVar1 = posConVar;
										}
										if (step == connectingPath.size() - 1 && linkedBlock == block2) {
											varInMBVar2 = linkedVar;
										} else if (step == connectingPath.size() - 1) {
											varInMBVar2 = posConVar;
										}

										ConfidencePair<Double, Integer> calculatedConf = calcConfCleanRunLinkedPreceding(
												linkedVar, posConVar, linkedBlock);
										resultConfidence.setConfidence(
												resultConfidence.confidence() * calculatedConf.confidence());
										resultConfidence.setSteps(resultConfidence.steps() + calculatedConf.steps());
									}
								}
							}
							block1 = block2;
						}
						if (varInMBVar1 == null || varInMBVar2 == null
								|| varInMBVar1.getAttributeAsInt("mergedBlocks", -1) != var1MB
								|| varInMBVar2.getAttributeAsInt("mergedBlocks", -1) != var2MB) {
							throw new Exception("Couldn't correctly find first and final connecting vars");
						}

						ConfidencePair<Double, Integer> tempConf;
						tempConf = calcConfCleanRunLinkedPreceding(var1, varInMBVar1, 0);
						resultConfidence.setConfidence(resultConfidence.confidence() * tempConf.confidence());
						resultConfidence.setSteps(resultConfidence.steps() + tempConf.steps());

						tempConf = calcConfCleanRunLinkedPreceding(var2, varInMBVar2, 0);
						resultConfidence.setConfidence(resultConfidence.confidence() * tempConf.confidence());
						resultConfidence.setSteps(resultConfidence.steps() + tempConf.steps());

						return resultConfidence;
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
					confP1 = calcConfCleanRunLinkedPreceding(var1, earliestTrioVar, 0);
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
					confP2 = calcConfCleanRunLinkedPreceding(var2, earliestTrioVar, 0);
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

	private ArrayList<Integer> calculateConnectingPath(int current, int destination,
			ArrayList<Integer[]> blockConnectionsList, ArrayList<Integer> path) {
		path.add(current);
		if (current == destination) {
			return path;
		}
		for (Integer[] connection : blockConnectionsList) {
			if (connection[0] == current && !path.contains(connection[1])) {
				int next = connection[1];
				ArrayList<Integer> result = calculateConnectingPath(next, destination, blockConnectionsList,
						new ArrayList<>(path));
				if (result != null) {
					return result;
				}
			}
			if (connection[1] == current && !path.contains(connection[0])) {
				int next = connection[0];
				ArrayList<Integer> result = calculateConnectingPath(next, destination, blockConnectionsList,
						new ArrayList<>(path));
				if (result != null) {
					return result;
				}
			}
		}

		return null;
	}

	public static boolean deepContains(List<Integer[]> list, Integer[] probe) {
		for (Integer[] element : list) {
			if (Arrays.deepEquals(element, probe)) {
				return true;
			}
		}
		return false;
	}

	private ConfidencePair<Double, Integer> calcConfCleanRunLinkedPreceding(VariantContext vc1, VariantContext vc2,
			int forcedJump) throws Exception {
		double product = 1.0;

		if (vc1.getStart() == vc2.getStart() && vc1.getReference().equals(vc2.getReference())
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

			VariantContext linkerV = null;
			;
			if (forcedJump > 0) {
				linkerV = calcComplicatedLinkerVar(linkerBlockVC, forcedJump);
			} else {
				linkerV = calcLinkerVar(linkerBlockVC);
			}

			ConfidencePair<Double, Integer> confP1 = calcConfCleanRunLinkedPreceding(linkerV, linkerBlockVC,
					forcedJump);
			ConfidencePair<Double, Integer> confP2 = calcConfCleanRunLinkedPreceding(otherBlockVC,
					getSimVC((VariantContext) linkerV.getAttribute("linkedPreceding")), forcedJump);

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

	private VariantContext calcComplicatedLinkerVar(VariantContext linkerBlockVC, int forcedJump) throws Exception {
		for (VariantContext v : this.getAllVariants()) {
			if (v.hasAttribute("linkedPreceding")
					&& v.getAttributeAsInt("mergedBlocks", -3) == linkerBlockVC.getAttributeAsInt("mergedBlocks", -4)
					&& getSimVC(((VariantContext) v.getAttribute("linkedPreceding", null)))
							.getAttributeAsInt("mergedBlocks", -1) == forcedJump) {
				return v;
			}
		}
		throw new Exception("NO LINKER VAR FOUND. Searching in block of: " + linkerBlockVC.toStringDecodeGenotypes());
	}

	private VariantContext calcLinkerVar(VariantContext linkerBlockVC) throws Exception {
		for (VariantContext v : this.getAllVariants()) {
			if (v.hasAttribute("linkedPreceding")
					&& v.getAttributeAsInt("mergedBlocks", -3) == linkerBlockVC.getAttributeAsInt("mergedBlocks", -4)
					&& getSimVC(((VariantContext) v.getAttribute("linkedPreceding", null))).getAttributeAsInt(
							"mergedBlocks", -1) == (linkerBlockVC.getAttributeAsInt("mergedBlocks", -2) - 1)) {
				return v;
			}
		}
		throw new Exception("NO LINKER VAR FOUND. Searching in block of: " + linkerBlockVC.toStringDecodeGenotypes());
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
			if (posVC.getStart() == trioVar.getStart() && posVC.getReference().equals(trioVar.getReference())
					&& posVC.getAlternateAllele(0).equals(trioVar
							.getAlternateAllele(0))) {
				if (trioVar.getAttributeAsBoolean("Innocuous", false)) {
					this.replaceVariant(new VariantContextBuilder(posVC)
							.genotypes(new GenotypeBuilder(posVC.getGenotype(PATIENT_ID)).phased(true)
									.attribute("TrioConfidence",
											trioVar.getGenotype(PATIENT_ID).getAnyAttribute("TrioConfidence"))
									.make())
							.attribute("Innocuous", true).make(), posVC, Strand.STRAND1);					
					return true;

				} else {
					this.replaceVariant(new VariantContextBuilder(posVC)
							.genotypes(new GenotypeBuilder(posVC.getGenotype(PATIENT_ID)).phased(true)
									.attribute("TrioConfidence",
											trioVar.getGenotype(PATIENT_ID).getAnyAttribute("TrioConfidence"))
									.make())
							.make(), posVC, Strand.STRAND1);
					return true;
				}
			}
		}
		for (VariantContext posVC : strand2) {
			if (posVC.getStart() == trioVar.getStart() && posVC.getReference().equals(trioVar.getReference())
					&& posVC.getAlternateAllele(0).equals(trioVar
							.getAlternateAllele(0))) {
				if (trioVar.getAttributeAsBoolean("Innocuous", false)) {
					this.replaceVariant(new VariantContextBuilder(posVC)
							.genotypes(new GenotypeBuilder(posVC.getGenotype(PATIENT_ID)).phased(true)
									.attribute("TrioConfidence",
											trioVar.getGenotype(PATIENT_ID).getAnyAttribute("TrioConfidence"))
									.make())
							.attribute("Innocuous", true).make(), posVC, Strand.STRAND2);					
					return true;
				} else {
					this.replaceVariant(new VariantContextBuilder(posVC)
							.genotypes(new GenotypeBuilder(posVC.getGenotype(PATIENT_ID)).phased(true)
									.attribute("TrioConfidence",
											trioVar.getGenotype(PATIENT_ID).getAnyAttribute("TrioConfidence"))
									.make())
							.make(), posVC, Strand.STRAND2);
					return true;
				}
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
				// Guarantee that tripHets are labeled as not phased to prevent problems when calculating confidence
				this.replaceVariant(new VariantContextBuilder(posVC).genotypes(
					new GenotypeBuilder(posVC.getGenotype(PATIENT_ID)).phased(false).make()).attribute("Innocuous", true).make(), posVC, Strand.STRAND1);				
				return true;
			}
		}
		for (VariantContext posVC : strand2) {
			if (posVC.getStart() == trioVar.getStart()) {
				this.replaceVariant(new VariantContextBuilder(posVC).genotypes(
						new GenotypeBuilder(posVC.getGenotype(PATIENT_ID)).phased(false).make()).attribute("Innocuous", true).make(), posVC,
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

	public void updateHMC() {
		if (possibleHMBC > highestMergedBlockCounter) {
			highestMergedBlockCounter = possibleHMBC;
		}
	}

}
