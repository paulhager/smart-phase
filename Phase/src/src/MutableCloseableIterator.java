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

import java.util.LinkedList;
import java.util.NoSuchElementException;
import java.util.Queue;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;

@SuppressWarnings("hiding")
public class MutableCloseableIterator<VariantContext> implements CloseableIterator<VariantContext> {
	
	private Queue<VariantContext> buffer = new LinkedList<>();
	
	public MutableCloseableIterator(CloseableIterator<VariantContext> closI) {
		buffer = new LinkedList<>(closI.toList());
	}
	
	public MutableCloseableIterator() {
		super();
	}
	
	
	@Override
	public VariantContext next() {
		if (!this.hasNext()) {
			throw new NoSuchElementException("Nothing left");
		}
		return this.buffer.poll();
	}

	@Override
	public void close() {}

	@Override
	public boolean hasNext() {
		if (this.buffer.isEmpty()){
			this.close();
			return false;
		}
		return true;
	}
	
	public void add(VariantContext vc) {
		buffer.add(vc);
	}

}
