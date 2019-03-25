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

import htsjdk.variant.variantcontext.VariantContext;
import smartPhase.HaplotypeBlock.Strand;

public class VariantStrandPair {

	private VariantContext outerVariant;
	private VariantContext innerVariant;
	private Strand innerVariantStrand;
	private Strand outerVariantStrand;
	
	public VariantStrandPair(VariantContext outerVariant, VariantContext innerVariant, Strand outerVariantStrand, Strand innerVariantStrand) {
		this.setOuterVariant(outerVariant);
		this.setInnerVariant(innerVariant);
		this.setOuterVariantStrand(outerVariantStrand);
		this.setInnerVariantStrand(innerVariantStrand);
	}

	public VariantContext getOuterVariant() {
		return outerVariant;
	}

	public void setOuterVariant(VariantContext outerVariant) {
		this.outerVariant = outerVariant;
	}

	public VariantContext getInnerVariant() {
		return innerVariant;
	}

	public void setInnerVariant(VariantContext innerVariant) {
		this.innerVariant = innerVariant;
	}

	public Strand getInnerVariantStrand() {
		return innerVariantStrand;
	}

	public void setInnerVariantStrand(Strand innerVariantStrand) {
		this.innerVariantStrand = innerVariantStrand;
	}

	public Strand getOuterVariantStrand() {
		return outerVariantStrand;
	}

	public void setOuterVariantStrand(Strand outerVariantStrand) {
		this.outerVariantStrand = outerVariantStrand;
	}
}
