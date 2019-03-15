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

import java.util.Iterator;
import java.util.LinkedList;
import java.util.NoSuchElementException;
import java.util.Queue;


@SuppressWarnings("hiding")
public class VariantIterator<SimpleVariant> implements Iterator<SimpleVariant>{
	private final Queue<SimpleVariant> buffer = new LinkedList<>();
	
	public VariantIterator(Iterator<SimpleVariant> iter) {
		super();
	}
	
	@Override
	public boolean hasNext() {
		if (this.buffer.isEmpty()){
			return false;
		}
		return true;
	}

	@Override
	public SimpleVariant next() {
		if (!this.hasNext()) {
		      throw new NoSuchElementException("Nothing left");
		    }
		    return this.buffer.poll();
	}
	
	public int currentSize(){
		return this.buffer.size();
	}
	
	

}
