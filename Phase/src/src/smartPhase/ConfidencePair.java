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

public class ConfidencePair<D, I> {
	
	private D confidence;
	private I steps;
	
	public ConfidencePair(D confidence, I steps){
		this.confidence = confidence;
		this.steps = steps;
	}
	
	public D confidence(){
		return this.confidence;
	}
	
	public I steps(){
		return this.steps;
	}
	
	public void setConfidence(D newConf){
		this.confidence = newConf;
	}
	
	public void setSteps(I newSteps){
		this.steps = newSteps;
	}
}
