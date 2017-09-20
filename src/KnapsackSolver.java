import java.util.Arrays;
import java.util.Comparator;

/**
 * An instance is a solver object for a specific instance of the "knapsack problem."
 * The knapsack problem is defined as follows:
 * 		Imagine you are robbing a house.
 * 		You have n items of varying value and weight that you can steal.
 * 		You have a knapsack that holds a certain maximum total weight.
 * 		Maximize the value you are carrying in the knapsack without exceeding 
 * 		the knapsack's weight capacity.
 * 
 * This solver class will provide an approximate solution to the given knapsack problem.
 * The approximation coefficient beta provides a lower bound on the total value of the 
 * items in the approximate solution. i.e. the approximate solution will have total value
 * of at least beta times the optimal solution, where beta is a double between 0 and 1.
 * 
 * The runtime efficiency is dependent on the square of the number of items, the max number of bits in
 * the binary representation of any of the numbers representing values or weights, and the coefficient beta.
 * In big-O notation, this algorithm takes time O(n^2 * s * beta), where n is the number of items, s is the 
 * max number of bits in the binary representation of any of the numbers, and beta is the approximation
 * coefficient.
 * 
 * To use this class, one must instantiate a KnapsackSolver object (thereby specifying all the 
 * necessary parameters of the knapsack problem). Then the user must run the solve() method on
 * that KnapsackSolver instance. After solving, the instance will contain the solution information,
 * namely the total value and total weight of the items in the approximate solution and a boolean[] 
 * specifying which items are in the approximate solution. Access these values using the methods
 * solvedValue(), solvedWeight() and solvedItems() respectively.
 * 
 * The user is also allowed to modify the capacity and beta parameters of an existing KnapsackSolver
 * instance. This is helpful if, for instance, a user wants to experiment with different capacities 
 * or beta values for an otherwise identical problem. After one of these parameters is changed, any
 * solution information previously stored in the KnapsackSolver instance will be erased since the 
 * problem has been changed and must be solved again.
 * 
 * @author Zach Mayle
 *
 */
public class KnapsackSolver {
	
	/**
	 * CLASS INVARIANTS
	 * 	
	 * 	1) the length of v and length of w must be equal
	 * 	2) n must equal the length of v (and therefore, length of w)
	 * 	3) solved is true iff the given knapsack problem instance is solved and this KnapsackSolver
	 * 		instance holds the solution information for the problem it represents
	 * 	4) if solved = false, then total_v = -1 and total_w = -1 and picked = {}
	 * 	5) if solved = true, then total_v = total value of the approximate solution
	 * 		and total_w = total weight of the approximate solution
	 * 		and picked = boolean array indicating which items are chosen
	 * 	
	 */
	
	private int n;  			// number of items to choose from
	private int[] v; 			// array holding value of each item ( v[i] = value of item i)
	private int[] w;			// array holding weight of each item ( w[] = weight of item i)
	private int W;				// weight capacity of the knapsack
	private double beta;		// approximation value
	private boolean[] picked;	// indicates which items are chosen in the solution, empty if no solution yet
	private int total_v;		// total value of the solution, -1 if no solution yet
	private int total_w;		// total weight of the solution, -1 if no solution yet
	private boolean solved;		// whether the knapsack problem has been solved for this system yet
	
	
	/** 
	 * Constructor: Creates a solver object for the specified instance of the Knapsack Problem.
	 * @param values: An array holding the value of each item to be considered for the knapsack. 
	 * 	Must have the same length as weights[].
	 * 	i.e. values[i] = value of item i. 
	 * @param weights: An array holding the weight of each item to be considered for the knapsack.
	 * 	Must have the same length as values[].
	 * 	i.e. weights[i] = weight of item i.
	 * @param capacity: The max total weight that the knapsack can hold. Must be > 0.
	 * @param beta: The approximation coefficient for this system. Must be between 0 and 1, exclusive.
	 * 	i.e. 0 < beta < 1
	 * 	The larger this beta value is, the longer the algorithm will take to compute.
	 */
	public KnapsackSolver(int[] values, int[] weights, int capacity, double beta) {
		if ( values.length != weights.length) {
			throw new IllegalArgumentException("The values list and weights list must have equal length.");
		}
		if (capacity <= 0) {
			throw new IllegalArgumentException("capacity must be an int > 0.");
		}
		if (beta <= 0 || beta >= 1) {
			throw new IllegalArgumentException("beta must be a double between 0 and 1. (0 < beta < 1)");
		}
		this.n = values.length;
		this.v = values;
		this.w = weights;
		this.W = capacity;
		this.beta = beta;
		this.picked = new boolean[this.n];
		this.total_v = -1;
		this.total_w = -1;
		this.solved = false;
	}
	
	
	// GETTER METHODS
	
	/**
	 * Getter for the number of items.
	 * @return Returns an int, the number of items considered in this instance of the Knapsack Problem.
	 */
	public int getNumberItems() {
		return n;
	}
	
	/**
	 * Getter for the item values.
	 * @return Returns an int[], representing the value of each item in the problem.
	 */
	public int[] getValues() {
		return v;
	}
	
	/**
	 * Getter for the item weights.
	 * @return Returns an int[], representing the weight of each item in the problem.
	 */
	public int[] getWeights() {
		return w;
	}
	
	/**
	 * Getter for the weight capacity.
	 * @return Returns an int, the max total weight the knapsack can hold in this problem..
	 */
	public int getCapacity() {
		return W;
	}
	
	/**
	 * Getter for the approximation coefficient.
	 * @return Returns a double, the approximation coefficient currently chosen for this problem.
	 */
	public double getBeta() {
		return beta;
	}
	
	
	// SETTER METHODS
	
	/**
	 * Setter for the weight capacity. On setting this value, the instance of the knapsack problem is
	 * modified and must be solved again using method solve() before you can access the new solution.
	 * @param capacity: int > 0. The max total weight capacity of the knapsack. 
	 */
	public void setCapacity(int capacity) {
		if (capacity <= 0) {
			throw new IllegalArgumentException("capacity must be an int > 0.");
		}
		this.W = capacity;
		this.solved = false;
		this.picked = new boolean[this.n];
		this.total_v = -1;
		this.total_w = -1;
	}
	
	/**
	 * Setter for the approximation coefficient beta. On setting this value, the instance of the 
	 * knapsack problem is modified and must be solved again using method solve() before you can
	 * access the new solution.
	 * @param capacity: double > 0 and < 1. The desired approximation coefficient for this problem.
	 */
	public void setBeta(double beta) {
		if (beta <= 0 || beta >= 1) {
			throw new IllegalArgumentException("beta must be a double between 0 and 1. (0 < beta < 1)");
		}
		this.beta = beta;
		this.solved = false;
		this.picked = new boolean[this.n];
		this.total_v = -1;
		this.total_w = -1;
	}
	
	
	// SOLUTION METHODS
	
	/**
	 * Returns the items chosen to be put in the knapsack in the current solution to this knapsack problem.
	 * @return boolean[] picked, where picked[i] is true iff item i is chosen for the given solution
	 * @exception IllegalStateException when this method is called before the given knapsack problem is solved.
	 */
	public boolean[] solvedItems() {
		if (!solved) {
			throw new IllegalStateException("The system has not yet been solved. "
					+ "Solve the system with the solve() method before attempting to access "
					+ "the items in the solution.");
		}
		return picked;
	}
	
	/**
	 * Returns the total value of all the items in the knapsack in the solution to this knapsack problem.
	 * @return int, the total value of items in the knapsack
	 * @exception IllegalStateException when this method is called before the given knapsack problem is solved.
	 */
	public int solvedValue() {
		if (!solved) {
			throw new IllegalStateException("The system has not yet been solved. "
					+ "Solve the system with the solve() method before attempting to access "
					+ "the total value of the solution.");
		}
		return total_v;
	}
	
	/**
	 * Returns the total weight of all the items in the knapsack in the solution to this knapsack problem.
	 * @return int, the total weight of items in the knapsack
	 * @exception IllegalStateException when this method is called before the given knapsack problem is solved.
	 */
	public int solvedWeight() {
		if (!solved) {
			throw new IllegalStateException("The system has not yet been solved. "
					+ "Solve the system with the solve() method before attempting to access "
					+ "the total weight of the solution.");
		}
		return total_w;
	}
	
	/**
	 * Solves the given instance of the knapsack problem approximately according to the given
	 * approximation coefficient beta and makes the solution available for access via the
	 * methods solvedItems(), solvedValue(), and solvedWeight().
	 * 
	 * If beta is the approximation coefficient for this KnapsackSolver instance, then 
	 * the approximate solution will have a total value of at least beta * value of optimal solution.
	 * 
	 */
	public void solve() {
		// EMGA (chooses most valuable items)
        Integer[] vSorted = sortedIndices(v);
        int weightEMGA = 0;
        int valueEMGA = 0;
        for (int i = vSorted.length -1; i >= 0; i--) {
        	if ( w[vSorted[i]] + weightEMGA <= W ) {
        		valueEMGA += v[vSorted[i]];
        		weightEMGA += w[vSorted[i]];
        	}
        	else break;
        }
        // now valueEMGA contains the value of the set formed by EMGA
        
        
        // GA (chooses items with best value to weight ratio)
        double[] rho = new double[n];
        for (int i = 0; i < n; i++) {
        	double vi = v[i];
        	double wi = w[i];
        	rho[i] = vi / wi;
        }
        int weightGA = 0;
        int valueGA = 0;
        Integer[] rhoSorted = sortedIndices(rho);
        for (int i = rhoSorted.length - 1; i >= 0; i--) {
        	if ( w[rhoSorted[i]] + weightGA <= W ) {
        		valueGA += v[rhoSorted[i]];
        		weightGA += w[rhoSorted[i]];
        	}
        	else break;
        }
        // now valueGA contains the value of the set formed by GA
        
        
        // find the max greedy value from the two greedy algorithms above
        int maxGreedy = Integer.max(valueGA, valueEMGA);
        
        
        // define D based on maxGreedy
        double epsilon = (1/beta) - 1;
        double D = epsilon * maxGreedy / (n * (epsilon + 1));
        
        
        // modify values
        int vt[] = new int[n];
        for (int i=0; i < n; i++) {
        	vt[i] = (int) Math.ceil(v[i]/D); // note how I defined vt differently from lecture
        }
        
        
        // define Vtilda
        int Vtilda = (int) Math.ceil((2 * n * (epsilon + 1))/ epsilon) + n;
        
        
        // initialize U, U[i,k] := min weight to fill the knapsack with value k using
        // only items from the set 1...i
        int U[][] = new int[n+1][Vtilda+1];
        U[0][0] = 0;
        for (int k = 1; k <= Vtilda; k++) {
        	U[0][k] = Integer.MAX_VALUE;
        }
        
        
        // Fill out U using dynamic programming.
        // Note that arrays v and w are zero indexed so I have to subtract 1 when searching
        // for the value or weight of the ith item.
        for (int i = 1; i <= n; i++) {
        	for (int k = 0; k <= Vtilda; k++) {
        		if ( k >= vt[i-1] ) {
        			int w1 = U[i-1][k];
        			int w2 = Integer.max(U[i-1][k-vt[i-1]] + w[i-1], U[i-1][k-vt[i-1]]);
        			if ( w1 <= w2 ) U[i][k] = w1;
        			else U[i][k] = w2;
        		} else {
        			U[i][k] = U[i-1][k];
        		}
        	}
        }
        
        
        // find max possible value kmax such that the weight constraint isn't violated
        int kmax = 0;
        for (int k = 0; k <= Vtilda; k++) {
        	if ( U[n][k] <= W ) {
        		if ( k > kmax ) kmax = k;
        	}
        }
        
        
        // trace back through U to construct the set that optimizes the rounded values
        int i = n;
        int k = kmax;
        while ( i > 0 && k > 0 ) {
        	if ( w[i-1] == U[i][k] ) {
        		if ( vt[i-1] >= k ) {
        			// then i is the only item I still need to put into S
        			picked[i-1] = true;
        		} else {
        			// then i is not in S
        			picked[i-1] = false;
        		}
        	} else if ( w[i-1] > U[i][k] ) {
        		// then i does not fit in S
        		picked[i-1] = false;
        	} else if ( U[i][k] < U[i-1][k] ) {
        		// then i has to be in S
        		picked[i-1] = true;
        	} else {
        		// i is not in S
        		picked[i-1] = false;
        	}
        	
        	// decrement k and i according to whether I put item i in the knapsack
        	// (I only decrement k if I put i in the knapsack)
        	if (picked[i-1]) k -= vt[i-1];
    		i -= 1;	
        }
        
        this.total_v = 0;
        this.total_w = 0;
        for (int j = 0; j < this.n; j++) {
        	if (picked[j]) {
        		this.total_v += v[j];
        		this.total_w += w[j];
        	}
        }
		
		this.solved = true;
	}
	
	
	// HELPER METHODS (given by my Algorithms professors Robert Kleinberg and Frans Schalekamp)
	private Integer[] sortedIndices( final double[] a )
    {
        Integer[] idx = new Integer[ n ];
        for ( int i = 0 ; i < a.length ; i++ )
        	idx[ i ] = new Integer( i );
        Arrays.sort
        ( idx, new Comparator<Integer>() 
	        {
			    @Override public int compare( final Integer i1, final Integer i2) 
			    {
			        return Double.compare( a[ i1 ], a[ i2 ] );
			    }
		    }
	    );
	    return idx;
    }
    
    private Integer[] sortedIndices( final int[] a )
    {
        Integer[] idx = new Integer[ n ];
        for ( int i = 0 ; i < a.length ; i++ )
        	idx[ i ] = new Integer( i );
        Arrays.sort
        ( idx, new Comparator<Integer>() 
	        {
			    @Override public int compare( final Integer i1, final Integer i2) 
			    {
			        return Integer.compare( a[ i1 ], a[ i2 ] );
			    }
		    }
	    );
	    return idx;
    }
}
