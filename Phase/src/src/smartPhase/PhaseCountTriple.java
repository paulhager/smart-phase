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

import java.util.Set;

import htsjdk.variant.variantcontext.VariantContext;


public class PhaseCountTriple <T, U> {
	private final Set<VariantContext> variants;
	private final SmartPhase.Phase phase;
	
	
	public PhaseCountTriple(Set<VariantContext> key, SmartPhase.Phase phase){
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
		PhaseCountTriple<Set<VariantContext>, SmartPhase.Phase> other = (PhaseCountTriple<Set<VariantContext>, SmartPhase.Phase>) obj;
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
