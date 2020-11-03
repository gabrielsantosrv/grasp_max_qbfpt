package problems.qbfpt;

import java.io.IOException;

/**
 * Class representing the inverse of the Quadractic Binary Function
 * ({@link QBFPT}), which is used since the GRASP is set by
 * default as a minimization procedure.
 * 
 * @author ccavellucci, fusberti
 */
public class QBFPT_Inverse extends QBFPT {

	/**
	 * Constructor for the QBFPT_Inverse class.
	 * 
	 * @param filename
	 *            Name of the file for which the objective function parameters
	 *            should be read.
	 * @throws IOException
	 *             Necessary for I/O operations.
	 */
	public QBFPT_Inverse(String filename) throws IOException {
		super(filename);
	}


	/* (non-Javadoc)
	 * @see problems.qbf.QBFPT#evaluate()
	 */
	@Override
	public Double evaluateQBFPT() {
		return -super.evaluateQBFPT();
	}
	
	/* (non-Javadoc)
	 * @see problems.qbf.QBFPT#evaluateInsertion(int)
	 */
	@Override
	public Double evaluateInsertionQBFPT(int i) {
		return -super.evaluateInsertionQBFPT(i);
	}
	
	/* (non-Javadoc)
	 * @see problems.qbf.QBFPT#evaluateRemoval(int)
	 */
	@Override
	public Double evaluateRemovalQBFPT(int i) {
		return -super.evaluateRemovalQBFPT(i);
	}
	
	/* (non-Javadoc)
	 * @see problems.qbf.QBFPT#evaluateExchange(int, int)
	 */
	@Override
	public Double evaluateExchangeQBFPT(int in, int out) {
		return -super.evaluateExchangeQBFPT(in,out);
	}

}
