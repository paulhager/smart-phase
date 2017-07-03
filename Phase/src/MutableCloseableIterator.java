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
