// Simple variant with chromosome string and start position
public class SimpleVariant {
	
	private String chrom;
	private int start;
	private int end;
	
	public SimpleVariant(String chrom, int start, int end){
		this.chrom = chrom;
		this.start = start;
		this.end = end;
	}
	
	
	
	public void setChrom(String chrom){
		this.chrom = chrom;
	}
	
	public String getChrom(){
		return this.chrom;
	}
	
	
	
	public void setStart(int start){
		this.start = start;
	}
	
	public int getStart(){
		return this.start;
	}
	
	
	
	public void setEnd(int end){
		this.end = end;
	}
	
	public int getEnd(){
		return this.end;
	}

}
