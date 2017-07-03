import java.util.Set;

import htsjdk.variant.variantcontext.VariantContext;


public class PhaseCountTriple <T, U> {
	private final Set<VariantContext> variants;
	private final smartPhase.Phase phase;
	
	
	public PhaseCountTriple(Set<VariantContext> key, smartPhase.Phase phase){
		this.phase = phase;
		this.variants = key;
	}
	
	/* (non-Javadoc)
	 * @see java.lang.Object#hashCode()
	 */
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((phase == null) ? 0 : phase.hashCode());
		result = prime * result + ((variants == null) ? 0 : variants.hashCode());
		return result;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (!(obj instanceof PhaseCountTriple))
			return false;
		@SuppressWarnings("unchecked")
		PhaseCountTriple<Set<VariantContext>, smartPhase.Phase> other = (PhaseCountTriple<Set<VariantContext>, smartPhase.Phase>) obj;
		if (phase != other.phase)
			return false;
		if (variants == null) {
			if (other.variants != null)
				return false;
		} else if (!variants.equals(other.variants))
			return false;
		return true;
	}
}
