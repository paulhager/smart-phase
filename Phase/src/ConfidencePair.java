
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
}
