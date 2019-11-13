package smartPhase;

import htsjdk.variant.variantcontext.VariantContext;

public class VariantTuple<V, T> {
	private final VariantContext v1;
	private final VariantContext v2;
	
	public VariantTuple(VariantContext v1, VariantContext v2) {
		this.v1 = v1;
		this.v2 = v2;
	}
	
	/* (non-Javadoc)
	 * @see java.lang.Object#hashCode()
	 */
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((v1 == null) ? 0 : v1.hashCode());
		result = prime * result + ((v2 == null) ? 0 : v2.hashCode());
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
		if (!(obj instanceof VariantTuple))
			return false;
		@SuppressWarnings("unchecked")
		VariantTuple<VariantContext, VariantContext> other = (VariantTuple<VariantContext, VariantContext>) obj;
		if (v1 != other.v1)
			return false;
		if (v2 != other.v2) {
			return false;
		} 
		return true;
	}

}
