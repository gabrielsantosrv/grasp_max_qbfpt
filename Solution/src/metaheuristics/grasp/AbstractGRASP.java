/**
 * 
 */
package metaheuristics.grasp;

import java.util.ArrayList;
import java.util.Random;

import problems.Evaluator;
import solutions.Solution;

/**
 * Abstract class for metaheuristic GRASP (Greedy Randomized Adaptive Search
 * Procedure). It consider a minimization problem.
 * 
 * @author ccavellucci, fusberti, vferrari, gabrielsantosrv, satoru27
 * @param <E>
 *            Generic type of the element which composes the solution.
 */
public abstract class AbstractGRASP<E> {

	public enum Construction {
		DEF,
		RPG
	}

	/**
	 * Value to represent the construction type.
	 * Can be default (DEF) or random plus greedy (RPG).
	 */
	private final Construction constructionType;

	/**
	 * Value that represents after how many iterations the random construction must turn greedy.
	 */
	private final int rpgP;

	/**
	 * flag that indicates whether the code should print more information on
	 * screen
	 */
	public static boolean verbose = true;

	/**
	 * a random number generator
	 */
	protected static Random rng = new Random(0);

	/**
	 * the objective function being optimized
	 */
	protected Evaluator<E> ObjFunction;

	/**
	 * the GRASP greediness-randomness parameter
	 */
	protected Double alpha;

	/**
	 * the best solution cost
	 */
	protected Double incumbentCost;

	/**
	 * the current solution cost
	 */
	protected Double currentCost;

	/**
	 * the best solution
	 */
	protected Solution<E> incumbentSol;

	/**
	 * the incumbent solution
	 */
	protected Solution<E> currentSol;

	/**
	 * the number of iterations the GRASP main loop executes.
	 */
	protected Integer iterations;

	/**
	 * the Candidate List of elements to enter the solution.
	 */
	protected ArrayList<E> CL;

	/**
	 * the Restricted Candidate List of elements to enter the solution.
	 */
	protected ArrayList<E> RCL;

	/**
	 * Creates the Candidate List, which is an ArrayList of candidate elements
	 * that can enter a solution.
	 * 
	 * @return The Candidate List.
	 */
	public abstract ArrayList<E> makeCL();

	/**
	 * Creates the Restricted Candidate List, which is an ArrayList of the best
	 * candidate elements that can enter a solution. The best candidates are
	 * defined through a quality threshold, delimited by the GRASP
	 * {@link #alpha} greedyness-randomness parameter.
	 * 
	 * @return The Restricted Candidate List.
	 */
	public abstract ArrayList<E> makeRCL();

	/**
	 * Updates the Candidate List according to the incumbent solution
	 * {@link #currentSol}. In other words, this method is responsible for
	 * updating which elements are still viable to take part into the solution.
	 */
	public abstract void updateCL();

	/**
	 * Creates a new solution which is empty, i.e., does not contain any
	 * element.
	 * 
	 * @return An empty solution.
	 */
	public abstract Solution<E> createEmptySol();

	/**
	 * The GRASP local search phase is responsible for repeatedly applying a
	 * neighborhood operation while the solution is getting improved, i.e.,
	 * until a local optimum is attained.
	 * 
	 * @return An local optimum solution.
	 */
	public abstract Solution<E> localSearch();
	
	/**
	 * Function to randomly choose a candidate from the RCL.
	 * 
	 * @return Random RCL candidate.
	 */
	public abstract E chooseRandom();

	/**
	 * Constructor for the AbstractGRASP class.
	 * 
	 * @param objFunction
	 *            The objective function being minimized.
	 * @param alpha
	 *            The GRASP greediness-randomness parameter (within the range
	 *            [0,1])
	 * @param iterations
	 *            The number of iterations which the GRASP will be executed.
	 */
	public AbstractGRASP(Evaluator<E> objFunction, Double alpha, Integer iterations, Construction constructionType, int rpgP) {
		this.ObjFunction = objFunction;
		this.alpha = alpha;
		this.iterations = iterations;
		this.constructionType = constructionType;
		this.rpgP = rpgP;
	}
	
	/**
	 * The GRASP constructive heuristic, which is responsible for building a
	 * feasible solution by selecting in a greedy-random fashion, candidate
	 * elements to enter the solution.
	 * 
	 * @return A feasible solution to the problem being minimized.
	 */
	public Solution<E> constructiveHeuristic() {
        int iter=0;

		CL = makeCL();
		RCL = makeRCL();
		currentSol = createEmptySol();
		currentCost = Double.POSITIVE_INFINITY;
		
		// Initialize alpha with random.
		if (constructionType == Construction.RPG){
            this.alpha = 1.0;
		}

		/* Main loop, which repeats until the stopping criteria is reached. */
		while (!constructiveStopCriteria()) {

			double maxCost = Double.NEGATIVE_INFINITY, minCost = Double.POSITIVE_INFINITY;
			currentCost = ObjFunction.evaluate(currentSol);
			updateCL();
			if(this.CL.size() == 0) break;

			/*
			 * Explore all candidate elements to enter the solution, saving the
			 * highest and lowest cost variation achieved by the candidates.
			 */
			for (E c : CL) {
				Double deltaCost = ObjFunction.evaluateInsertionCost(c, currentSol);
				if (deltaCost < minCost)
					minCost = deltaCost;
				if (deltaCost > maxCost)
					maxCost = deltaCost;
			}
            
            /* Random plus greedy.
             *  Iterations [0,p) - Random: alpha=1
             *	Iterations [p,n) - Greedy: alpha=0
             */
			if (constructionType == Construction.RPG && iter == this.rpgP){
                this.alpha = 0.0;
			}

			/*
			 * Among all candidates, insert into the RCL those with the highest
			 * performance using parameter alpha as threshold.
			 */
			for (E c : CL) {
				Double deltaCost = ObjFunction.evaluateInsertionCost(c, currentSol);
				if (deltaCost <= minCost + alpha * (maxCost - minCost)) {
					RCL.add(c);
				}
			}

			/* Choose a candidate randomly from the RCL */
			E inCand = chooseRandom();
			CL.remove(inCand);
			currentSol.add(inCand);
			ObjFunction.evaluate(currentSol);
			RCL.clear();
            
            // Increase iteration count.
            iter++;
		}

		return currentSol;
	}

	/**
	 * The GRASP mainframe. It consists of a loop, in which each iteration goes
	 * through the constructive heuristic and local search. The best solution is
	 * returned as result.
	 * 
	 * @return The best feasible solution obtained throughout all iterations.
	 */
	public Solution<E> solve(double maxTime) {
		int i;
		long startTime = System.currentTimeMillis();
		long endTime;
		double totalTime;
		incumbentSol = createEmptySol();
		
		for (i = 0; i < iterations; i++) {
			constructiveHeuristic();
			localSearch();
			if (incumbentSol.cost > currentSol.cost) {
				incumbentSol = new Solution<E>(currentSol);
				if (verbose)
					System.out.println("(Iter. " + i + ") BestSol = " + incumbentSol);
			}
			
			endTime   = System.currentTimeMillis();
			totalTime = (endTime - startTime)/(double)1000;

			//if it exceeded the time limit of 1800s (30 min), then break the loop
			if(totalTime > maxTime) break;
		}
		
		if(verbose)
			System.out.println("Total iterations: " + i);

		return incumbentSol;
	}

	/**
	 * A standard stopping criteria for the constructive heuristic is to repeat
	 * until the incumbent solution improves by inserting a new candidate
	 * element.
	 * 
	 * @return true if the criteria is met.
	 */
	public Boolean constructiveStopCriteria() {
		return (currentCost > currentSol.cost) ? false : true;
	}

}
