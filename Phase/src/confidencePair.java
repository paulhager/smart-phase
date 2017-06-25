
public class confidencePair<D, I> {
	
	private D confidence;
	private I steps;
	
	public confidencePair(D confidence, I steps){
		this.confidence = confidence;
		this.steps = steps;
	}
	
	public D confidence(){
		return this.confidence;
	}
	
	public I steps(){
		return this.steps;
	}
}
