package problems.qbfpt;

import problems.Evaluator;
import solutions.Solution;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * A quadractic binary function (QBFPT) is a function that can be expressed as the
 * sum of quadractic terms: f(x) = \sum{i,j}{a_{ij}*x_i*x_j}. In matrix form
 * a QBFPT can be expressed as f(x) = x'.A.x
 * The problem of minimizing a QBFPT is NP-hard [1], even when no constraints
 * are considered.
 * 
 * [1] Kochenberger, et al. The unconstrained binary quadratic programming
 * problem: a survey. J Comb Optim (2014) 28:58–81. DOI
 * 10.1007/s10878-014-9734-0.
 * 
 * @author ccavellucci, fusberti
 *
 */
public class QBFPT implements Evaluator<Integer> {

	private static final int GPI1 = 131;
	private static final int GPI2 = 1031;
	private static final int HPI1 = 193;
	private static final int HPI2 = 1093;

	/**
	 * Dimension of the domain.
	 */
	public final Integer size;

	/**
	 * The array of numbers representing the domain.
	 */
	public final Double[] variables;

	/**
	 * The matrix A of coefficients for the QBFPT f(x) = x'.A.x
	 */
	public Double[][] A;

	private ArrayList<int[]> triples;

	/**
	 * The constructor for QuadracticBinaryFunction class. The filename of the
	 * input for setting matrix of coefficients A of the QBFPT. The dimension of
	 * the array of variables x is returned from the {@link #readInput} method.
	 * 
	 * @param filename
	 *            Name of the file containing the input for setting the QBFPT.
	 * @throws IOException
	 *             Necessary for I/O operations.
	 */
	public QBFPT(String filename) throws IOException {
		size = readInput(filename);
		variables = allocateVariables();
		triples = generate_triples();
	}
	
	private int[] generate_triple_aux(int u, int n){
		int l,g,h;
		int aux;
		
		// Generating first l(u).
		l = 1 + ((GPI1 * (u - 1) + GPI2) % n);
		
		// Choosing g(u).
		g = (l != u) ? l:(1 + (l % n));

		// Generating second l(u).
		l = 1 + ((HPI1 * (u - 1) + HPI2) % n);
		
		// Choosing h(u).
		if(l != u && l != g){
			h = l;
		}else{
			aux = (1 + (l % n));
			if(aux != u && aux != g){
				h = aux;
			}else{
				h = 1 + ((l + 1) % n);
			}
		}
		
		// Building output
		int triple[] = {u-1, g-1, h-1};
		Arrays.sort(triple);
		return triple;
	}

	/**
	 * Generate forbidden triples for QBFPT.
	 * Triples are in the format sort{u, g(u), h(u)}
	 * Both functions are based on a l(u) function.
	 * l(u) = 1 + ((pi1 * (u-1) + pi2) mod n)
	 * g(u) = 1+l(u) mod n if l(u) = u, else l(u)
	 * h(u) = {l(u) if l(u) != u and l(u) != g(u)
	 * 		   1 + l(u) mod n if this value != u and g(u)
	 * 		   else 1 + (l(u)+1) mod n}
	 * Values of pi are different for g(u) and h(u). 
	 * @return _triples ArrayList of forbidden triples.
	 */
	private ArrayList<int[]> generate_triples(){
		ArrayList<int[]> _triples = new ArrayList<>();
		for(int i=1; i <= this.size; i++){
			_triples.add(generate_triple_aux(i, this.size));
		}
		return _triples;
	}

	/**
	 * Test whether the solution obtained with the insertion of i will be feasible.
	 * For this, checks triple restriction.
	 * Checks if the 2 other elements of all triples with i are in the current solution.
	 * @param i Value to check feasibility.
	 * @return feasible: true if i is feasible, false otherwise.
	 */
	public boolean is_feasible(int i){
		boolean feasible = true;
		double sum;
		
		// Check all triples.
		for(int [] tuple : this.triples){
			
			// If the value is in a triple, check if triple is in the current solution.
			if(i == tuple[0] || i == tuple[1] || i == tuple[2]){
				sum = 0;
				for(int k : tuple){
					sum += this.variables[k];
				}

				//if there are 2 indices in a tuple that have already been set to 1
				//and the index i isn't of theses indices, then set it to 1 will
				//complete a prohibited tuple.
				if(sum == 2.0 && this.variables[i] == 0.0){
					feasible = false;
					break;
				}
			}
		}

		return feasible;
	}

	/**
	 * Evaluates the value of a solution by transforming it into a vector. This
	 * is required to perform the matrix multiplication which defines a QBFPT.
	 * 
	 * @param sol
	 *            the solution which will be evaluated.
	 */
	public void setVariables(Solution<Integer> sol) {

		resetVariables();
		if (!sol.isEmpty()) {
			for (Integer elem : sol) {
				variables[elem] = 1.0;
			}
		}

	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see problems.Evaluator#getDomainSize()
	 */
	@Override
	public Integer getDomainSize() {
		return size;
	}

	/**
	 * {@inheritDoc} In the case of a QBFPT, the evaluation correspond to
	 * computing a matrix multiplication x'.A.x. A better way to evaluate this
	 * function when at most two variables are modified is given by methods
	 * {@link #evaluateInsertionQBFPT(int)}, {@link #evaluateRemovalQBFPT(int)} and
	 * {@link #evaluateExchangeQBFPT(int,int)}.
	 * 
	 * @return The evaluation of the QBFPT.
	 */
	@Override
	public Double evaluate(Solution<Integer> sol) {

		setVariables(sol);
		return sol.cost = evaluateQBFPT();

	}

	/**
	 * Evaluates a QBFPT by calculating the matrix multiplication that defines the
	 * QBFPT: f(x) = x'.A.x .
	 * 
	 * @return The value of the QBFPT.
	 */
	public Double evaluateQBFPT() {

		Double aux = (double) 0, sum = (double) 0;
		Double vecAux[] = new Double[size];

		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				aux += variables[j] * A[i][j];
			}
			vecAux[i] = aux;
			sum += aux * variables[i];
			aux = (double) 0;
		}

		return sum;

	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see problems.Evaluator#evaluateInsertionCost(java.lang.Object,
	 * solutions.Solution)
	 */
	@Override
	public Double evaluateInsertionCost(Integer elem, Solution<Integer> sol) {

		setVariables(sol);
		return evaluateInsertionQBFPT(elem);

	}

	/**
	 * Determines the contribution to the QBFPT objective function from the
	 * insertion of an element.
	 * 
	 * @param i
	 *            Index of the element being inserted into the solution.
	 * @return The variation of the objective function resulting from the
	 *         insertion.
	 */
	public Double evaluateInsertionQBFPT(int i) {

		if (variables[i] == 1)
			return 0.0;

		return evaluateContributionQBFPT(i);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see problems.Evaluator#evaluateRemovalCost(java.lang.Object,
	 * solutions.Solution)
	 */
	@Override
	public Double evaluateRemovalCost(Integer elem, Solution<Integer> sol) {

		setVariables(sol);
		return evaluateRemovalQBFPT(elem);

	}

	/**
	 * Determines the contribution to the QBFPT objective function from the
	 * removal of an element.
	 * 
	 * @param i
	 *            Index of the element being removed from the solution.
	 * @return The variation of the objective function resulting from the
	 *         removal.
	 */
	public Double evaluateRemovalQBFPT(int i) {

		if (variables[i] == 0)
			return 0.0;

		return -evaluateContributionQBFPT(i);

	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see problems.Evaluator#evaluateExchangeCost(java.lang.Object,
	 * java.lang.Object, solutions.Solution)
	 */
	@Override
	public Double evaluateExchangeCost(Integer elemIn, Integer elemOut, Solution<Integer> sol) {

		setVariables(sol);
		return evaluateExchangeQBFPT(elemIn, elemOut);

	}

	/**
	 * Determines the contribution to the QBFPT objective function from the
	 * exchange of two elements one belonging to the solution and the other not.
	 * 
	 * @param in
	 *            The index of the element that is considered entering the
	 *            solution.
	 * @param out
	 *            The index of the element that is considered exiting the
	 *            solution.
	 * @return The variation of the objective function resulting from the
	 *         exchange.
	 */
	public Double evaluateExchangeQBFPT(int in, int out) {

		Double sum = 0.0;

		if (in == out)
			return 0.0;
		if (variables[in] == 1)
			return evaluateRemovalQBFPT(out);
		if (variables[out] == 0)
			return evaluateInsertionQBFPT(in);

		sum += evaluateContributionQBFPT(in);
		sum -= evaluateContributionQBFPT(out);
		sum -= (A[in][out] + A[out][in]);

		return sum;
	}

	/**
	 * Determines the contribution to the QBFPT objective function from the
	 * insertion of an element. This method is faster than evaluating the whole
	 * solution, since it uses the fact that only one line and one column from
	 * matrix A needs to be evaluated when inserting a new element into the
	 * solution. This method is different from {@link #evaluateInsertionQBFPT(int)},
	 * since it disregards the fact that the element might already be in the
	 * solution.
	 * 
	 * @param i
	 *            index of the element being inserted into the solution.
	 * @return the variation of the objective function resulting from the
	 *         insertion.
	 */
	private Double evaluateContributionQBFPT(int i) {

		Double sum = 0.0;

		for (int j = 0; j < size; j++) {
			if (i != j)
				sum += variables[j] * (A[i][j] + A[j][i]);
		}
		sum += A[i][i];

		return sum;
	}

	/**
	 * Responsible for setting the QBFPT function parameters by reading the
	 * necessary input from an external file. This method reads the domain's
	 * dimension and matrix {@link #A}.
	 * 
	 * @param filename
	 *            Name of the file containing the input for setting the black
	 *            box function.
	 * @return The dimension of the domain.
	 * @throws IOException
	 *             Necessary for I/O operations.
	 */
	protected Integer readInput(String filename) throws IOException {

		Reader fileInst = new BufferedReader(new FileReader(filename));
		StreamTokenizer stok = new StreamTokenizer(fileInst);

		stok.nextToken();
		Integer _size = (int) stok.nval;
		A = new Double[_size][_size];

		for (int i = 0; i < _size; i++) {
			for (int j = i; j < _size; j++) {
				stok.nextToken();
				A[i][j] = stok.nval;
				if (j>i)
					A[j][i] = 0.0;
			}
		}

		return _size;

	}

	/**
	 * Reserving the required memory for storing the values of the domain
	 * variables.
	 * 
	 * @return a pointer to the array of domain variables.
	 */
	protected Double[] allocateVariables() {
		Double[] _variables = new Double[size];
		return _variables;
	}

	/**
	 * Reset the domain variables to their default values.
	 */
	public void resetVariables() {
		Arrays.fill(variables, 0.0);
	}

	/**
	 * Prints matrix {@link #A}.
	 */
	public void printMatrix() {

		for (int i = 0; i < size; i++) {
			for (int j = i; j < size; j++) {
				System.out.print(A[i][j] + " ");
			}
			System.out.println();
		}

	}

	public void printTriples(){
		for(int[] tuple : triples)
			System.out.print("("+(tuple[0]+1)+","+(tuple[1]+1)+","+(tuple[2]+1)+"), ");

		System.out.println();
	}
}
