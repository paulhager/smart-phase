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
