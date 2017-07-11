import java.util.ArrayList;
import java.util.HashMap;

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
		vc = new VariantContextBuilder(vc).attribute("mergedBlocks", 1).make();
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
	 * 			  - Variant to be replaced
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
	 */
	public void addVariantsMerge(ArrayList<VariantContext> vcList, Strand s, int mergeBlockCntr) {

		for (VariantContext vc : vcList) {
			
			vc = new VariantContextBuilder(vc).attribute("mergedBlocks", mergeBlockCntr).make();
			strandVariants.get(s).add(vc);

			if (vc.getStart() < blockStart) {
				blockStart = vc.getStart();
			}

			if (vc.getEnd() > blockEnd) {
				blockEnd = vc.getEnd();
			}
		}
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
			if (posVC.getStart() == vc.getStart() && posVC.getContig().equals(vc.getContig()) && posVC.getReference().equals(vc.getReference()) && posVC.getAlternateAllele(0).equals(vc.getAlternateAllele(0))) {
				return Strand.STRAND1;
			}
		}
		for (VariantContext posVC : strand2) {
			if (posVC.getStart() == vc.getStart() && posVC.getContig().equals(vc.getContig()) && posVC.getReference().equals(vc.getReference()) && posVC.getAlternateAllele(0).equals(vc.getAlternateAllele(0))) {
				return Strand.STRAND2;
			}
		}
		return null;
	}
	
	/**
	 * Grabs saved variantContext with just start position
	 * @param vc - VariantContext with same start position as the one in block to be grabbed
	 * @return VariantContext
	 */
	public VariantContext getSimVC(VariantContext vc){
		for (VariantContext posVC : strand1) {
			if (posVC.getStart() == vc.getStart() && posVC.getContig().equals(vc.getContig()) && posVC.getReference().equals(vc.getReference()) && posVC.getAlternateAllele(0).equals(vc.getAlternateAllele(0))) {
				return posVC;
			}
		}
		for (VariantContext posVC : strand2) {
			if (posVC.getStart() == vc.getStart() && posVC.getContig().equals(vc.getContig()) && posVC.getReference().equals(vc.getReference()) && posVC.getAlternateAllele(0).equals(vc.getAlternateAllele(0))) {
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
		confidencePair<Double, Integer> cpTrio1 = multiplyConfidence(vc1, nearestTrioVC1);
		cpTrio1.setConfidence(cpTrio1.confidence()*(double)nearestTrioVC1.getGenotype(PATIENT_ID).getAnyAttribute("TrioConfidence"));
		
		confidencePair<Double, Integer> cpTrio2 = multiplyConfidence(vc2, nearestTrioVC2);
		cpTrio1.setConfidence(cpTrio2.confidence()*(double)nearestTrioVC2.getGenotype(PATIENT_ID).getAnyAttribute("TrioConfidence"));
		
		// Try to calculate confidence using only read conf. scores if both variants are in same mergeBlock. 
		if(vc1.getAttributeAsInt("mergedBlocks", -1) == vc2.getAttributeAsInt("mergedBlocks", -1)){
			confidencePair<Double, Integer> cpRead = multiplyConfidence(vc1, vc2);
			
			// Compare amount of steps to determine which conf. score is returned
			if(cpRead.steps() < (cpTrio1.steps() + cpTrio2.steps())/2){
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

	private confidencePair<Double, Integer> multiplyConfidence(VariantContext vc1, VariantContext vc2) throws Exception {
		double product = 1;
		
		if(vc1 == vc2){
			return new confidencePair<Double, Integer>(1.0, 0);
		}
		
		// Variants must be in same mergeBlock, or else something is wrong
		if (vc1.getAttributeAsInt("mergedBlocks", -1) != vc2.getAttributeAsInt("mergedBlocks", -1)) {
			System.out.println(vc1.toStringDecodeGenotypes());
			System.out.println(vc2.toStringDecodeGenotypes());
			throw new Exception("No trio and not in same block.... Hmmm..... Something went wrong.");
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
		int cnt = 1;
		while (curVC.getStart() != end) {
			product = product * (double) curVC.getGenotype(PATIENT_ID).getAnyAttribute("ReadConfidence");
			cnt++;
			curVC = this.getSimVC((VariantContext) curVC.getAttribute("Preceding"));
		}
		product = product * (double) curVC.getGenotype(PATIENT_ID).getAnyAttribute("ReadConfidence");
		
				
		return new confidencePair<Double, Integer>(product, cnt);
	}

	/**
	 * Set phased attribute of a variant's patient_GT context to true.
	 * @param trioVar - Variant whose phase is to be changed
	 * @return True - If var found
	 * 		   False - If var not found
	 */
	public boolean setPhased(VariantContext trioVar) {
		for (VariantContext posVC : strand1) {
			if (posVC.getStart() == trioVar.getStart()) {
				if(trioVar.getAttributeAsBoolean("Innocuous", false)){
					this.replaceVariant(new VariantContextBuilder(posVC).genotypes(new GenotypeBuilder(posVC.getGenotype(PATIENT_ID)).phased(true).attribute("TrioConfidence", trioVar.getGenotype(PATIENT_ID).getAnyAttribute("TrioConfidence")).make()).attribute("Innocuous", true).make(), posVC, Strand.STRAND1);
				} else {
					this.replaceVariant(new VariantContextBuilder(posVC).genotypes(new GenotypeBuilder(posVC.getGenotype(PATIENT_ID)).phased(true).attribute("TrioConfidence", trioVar.getGenotype(PATIENT_ID).getAnyAttribute("TrioConfidence")).make()).make(), posVC, Strand.STRAND1);
				}
				return true;
			}
		}
		for (VariantContext posVC : strand2) {
			if (posVC.getStart() == trioVar.getStart()) {
				if(trioVar.getAttributeAsBoolean("Innocuous", false)){
					this.replaceVariant(new VariantContextBuilder(posVC).genotypes(new GenotypeBuilder(posVC.getGenotype(PATIENT_ID)).phased(true).attribute("TrioConfidence", trioVar.getGenotype(PATIENT_ID).getAnyAttribute("TrioConfidence")).make()).attribute("Innocuous", true).make(), posVC, Strand.STRAND2);
				} else {
					this.replaceVariant(new VariantContextBuilder(posVC).genotypes(new GenotypeBuilder(posVC.getGenotype(PATIENT_ID)).phased(true).attribute("TrioConfidence", trioVar.getGenotype(PATIENT_ID).getAnyAttribute("TrioConfidence")).make()).make(), posVC, Strand.STRAND2);
				}
				return true;
			}
		}
		return false;
	}

	
	/**
	 * Set Innocuous attribute of a variant's patient_GT context to true.
	 * @param trioVar - Variant whose phase is to be changed
	 * @return True - If var found
	 * 		   False - If var not found
	 */
	public boolean setTripHet(VariantContext trioVar) {
		for (VariantContext posVC : strand1) {
			if (posVC.getStart() == trioVar.getStart()) {
				this.replaceVariant(new VariantContextBuilder(posVC).attribute("Innocuous", true).make(), posVC, Strand.STRAND1);
				return true;
			}
		}
		for (VariantContext posVC : strand2) {
			if (posVC.getStart() == trioVar.getStart()) {
				this.replaceVariant(new VariantContextBuilder(posVC).attribute("Innocuous", true).make(), posVC, Strand.STRAND2);
				return true;
			}
		}
		return false;
	}

}
