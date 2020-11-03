package problems.qbfpt.solvers;

import metaheuristics.grasp.AbstractGRASP;
import problems.qbfpt.QBFPT;
import problems.qbfpt.QBFPT_Inverse;
import solutions.Solution;

import java.io.IOException;
import java.util.ArrayList;


/**
 * Metaheuristic GRASP (Greedy Randomized Adaptive Search Procedure) for
 * obtaining an optimal solution to a QBFPT (Quadractive Binary Function --
 * {@link #QuadracticBinaryFunction}). Since by default this GRASP considers
 * minimization problems, an inverse QBFPT function is adopted.
 * 
 * @author ccavellucci, fusberti
 */
public class GRASP_QBFPT extends AbstractGRASP<Integer> {
	private enum SearchStrategy {
		FI,
		BI
	}
	
	/**
	 * Value to represent local search type.
	 * Can be first-improving (FI) or best-improving (BI). 
	 */
	private final SearchStrategy searchType;
	
	/**
	 * Constructor for the GRASP_QBFPT class. An inverse QBFPT objective function is
	 * passed as argument for the superclass constructor.
	 * 
	 * @param alpha
	 *            The GRASP greediness-randomness parameter (within the range
	 *            [0,1])
	 * @param iterations
	 *            The number of iterations which the GRASP will be executed.
	 * @param filename
	 *            Name of the file for which the objective function parameters
	 *            should be read.
	 * @param type
	 * 			  Type of local search strategy to be used.
	 * @throws IOException
	 *             necessary for I/O operations.
	 */
	public GRASP_QBFPT(Double alpha, Integer iterations, String filename, SearchStrategy type) throws IOException {
		super(new QBFPT_Inverse(filename), alpha, iterations);
		this.searchType = type;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see grasp.abstracts.AbstractGRASP#makeCL()
	 */
	@Override
	public ArrayList<Integer> makeCL() {

		ArrayList<Integer> _CL = new ArrayList<Integer>();
		for (int i = 0; i < ObjFunction.getDomainSize(); i++) {
			Integer cand = i;
			_CL.add(cand);
		}

		return _CL;

	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see grasp.abstracts.AbstractGRASP#makeRCL()
	 */
	@Override
	public ArrayList<Integer> makeRCL() {
		return new ArrayList<Integer>();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see grasp.abstracts.AbstractGRASP#updateCL()
	 */
	@Override
	public void updateCL() {
		ArrayList<Integer> CL_copy = new ArrayList<>(this.CL);
		for (Integer c : CL_copy){
			if(!((QBFPT)this.ObjFunction).is_feasible(c)){
				this.CL.remove(c);
			}
		}
	}

	/**
	 * {@inheritDoc}
	 * 
	 * This createEmptySol instantiates an empty solution and it attributes a
	 * zero cost, since it is known that a QBFPT solution with all variables set
	 * to zero has also zero cost.
	 */
	@Override
	public Solution<Integer> createEmptySol() {
		Solution<Integer> sol = new Solution<Integer>();
		sol.cost = 0.0;
		return sol;
	}

	/**
	 * {@inheritDoc}
	 * 
	 * The local search operator developed for the QBFPT objective function is
	 * composed by the neighborhood moves Insertion, Removal and 2-Exchange.
	 */
	@Override
	public Solution<Integer> localSearch(){
		Solution<Integer> sol;
		
		// Check local search method.
		if (this.searchType == SearchStrategy.BI)
			sol = localSearchBestImproving();
		else
			sol = localSearchFirstImproving();
	
		return sol;
	}
	
	/**
	 * Best-Improving local search.
	 */
	private Solution<Integer> localSearchBestImproving() {

		Double minDeltaCost;
		Integer bestCandIn = null, bestCandOut = null;

		do {
			minDeltaCost = Double.POSITIVE_INFINITY;
			updateCL();
				
			// Evaluate insertions
			for (Integer candIn : CL) {
				double deltaCost = ObjFunction.evaluateInsertionCost(candIn, currentSol);
				if (deltaCost < minDeltaCost) {
					minDeltaCost = deltaCost;
					bestCandIn = candIn;
					bestCandOut = null;
				}
			}
			// Evaluate removals
			for (Integer candOut : currentSol) {
				double deltaCost = ObjFunction.evaluateRemovalCost(candOut, currentSol);
				if (deltaCost < minDeltaCost) {
					minDeltaCost = deltaCost;
					bestCandIn = null;
					bestCandOut = candOut;
				}
			}
			// Evaluate exchanges
			for (Integer candIn : CL) {
				for (Integer candOut : currentSol) {
					double deltaCost = ObjFunction.evaluateExchangeCost(candIn, candOut, currentSol);
					if (deltaCost < minDeltaCost) {
						minDeltaCost = deltaCost;
						bestCandIn = candIn;
						bestCandOut = candOut;
					}
				}
			}
			// Implement the best move, if it reduces the solution cost.
			if (minDeltaCost < -Double.MIN_VALUE) {
				if (bestCandOut != null) {
					currentSol.remove(bestCandOut);
					CL.add(bestCandOut);
				}
				if (bestCandIn != null) {
					currentSol.add(bestCandIn);
					CL.remove(bestCandIn);
				}
				ObjFunction.evaluate(currentSol);
			}
		} while (minDeltaCost < -Double.MIN_VALUE);

		return null;
	}
	
	/**
	 * First-Improving local search.
	 */
	private Solution<Integer> localSearchFirstImproving() {

		Double minDeltaCost;
		Integer firstCandIn = null, firstCandOut = null;

		do {
			minDeltaCost = Double.POSITIVE_INFINITY;
			updateCL();
				
			// Evaluate insertions
			for (Integer candIn : CL) {
				double deltaCost = ObjFunction.evaluateInsertionCost(candIn, currentSol);
				if (deltaCost < minDeltaCost) {
					minDeltaCost = deltaCost;
					firstCandIn = candIn;
					firstCandOut = null;
					break;
				}
			}
			
			// Evaluate removals
			for (Integer candOut : currentSol) {
				double deltaCost = ObjFunction.evaluateRemovalCost(candOut, currentSol);
				if (deltaCost < minDeltaCost) {
					minDeltaCost = deltaCost;
					firstCandIn = null;
					firstCandOut = candOut;
					break;
				}
			}
			
			// Evaluate exchanges
			for (Integer candIn : CL) {
				for (Integer candOut : currentSol) {
					double deltaCost = ObjFunction.evaluateExchangeCost(candIn, candOut, currentSol);
					if (deltaCost < minDeltaCost) {
						minDeltaCost = deltaCost;
						firstCandIn = candIn;
						firstCandOut = candOut;
						break;
					}
				}
			}
			// Implement the best first move, if it reduces the solution cost.
			if (minDeltaCost < -Double.MIN_VALUE) {
				if (firstCandOut != null) {
					currentSol.remove(firstCandOut);
					CL.add(firstCandOut);
				}
				if (firstCandIn != null) {
					currentSol.add(firstCandIn);
					CL.remove(firstCandIn);
				}
				ObjFunction.evaluate(currentSol);
			}
		} while (minDeltaCost < -Double.MIN_VALUE);

		return null;
	}
	
	/**
	 * {@inheritDoc}
	 * 
	 * The QBFPT random choice follows a bias function.
	 * Get the bias for each value, calculate probability, and choose element.
	 */
	@Override
	public Integer chooseRandom(){
		double probs[] = new double[RCL.size()];
		double totalBias = 0;
		int i;
		
		// Get bias.
		for(i=0; i<RCL.size(); i++) {
			probs[i] = 1;
			totalBias += 1;
		}
		
		// Calculate probabilities.
		for(i=0; i<RCL.size(); i++) {
			probs[i] /= totalBias;
		}
		
		// Get random value from weighted probs.
		int rndIndex=0;
		double rndValue = rng.nextDouble();
		for(i=0; i<RCL.size(); i++) {
			if(rndValue < probs[i]) {
				rndIndex = i;
				break;
			}
			rndValue -= probs[i];
		}
		
		return RCL.get(rndIndex);
	}

	/**
	 * A main method used for testing the GRASP metaheuristic.
	 */
	public static void main(String[] args) throws IOException {

		long startTime = System.currentTimeMillis();
		GRASP_QBFPT grasp = new GRASP_QBFPT(0.25, 
											1000,
											"instances/qbf020", 
											SearchStrategy.BI);
		
		Solution<Integer> bestSol = grasp.solve();
		System.out.println("maxVal = " + bestSol);

		long endTime   = System.currentTimeMillis();
		long totalTime = endTime - startTime;
		System.out.println("Time = "+(double)totalTime/(double)1000+" seg");
	}
}
