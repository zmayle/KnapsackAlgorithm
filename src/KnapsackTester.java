import static org.junit.Assert.*;

import org.junit.Test;

public class KnapsackTester {

	@Test
	public void test() {
		int[] values = {0, 1, 3, 5, 10, 9, 8, 7, 6, 3, 2, 8};
		int[] weight = {1, 4, 5, 2, 12, 1, 7, 6, 2, 8, 12, 8};
		int cap = 18;
		double beta = 0.99;
		KnapsackSolver knapsack = new KnapsackSolver(values, weight, cap, beta);
		knapsack.solve();
		int v = knapsack.solvedValue();
		int w = knapsack.solvedWeight();
		boolean[] picked = knapsack.solvedItems();
		System.out.println("value: " + v);
		System.out.println("weight: " + w);
		
		int a = 0;
		for (int i = 0; i < picked.length; i++) {
			if (picked[i]) {
				if (a > 0) {
				   System.out.print(", ");
				}
				System.out.print(values[i]);
				a += 1;
			}
		}
		
		knapsack.setBeta(0.9999);
		knapsack.solve();
		System.out.println("value: " + knapsack.solvedValue());
		System.out.println("weight: " + knapsack.solvedWeight());
	}

}
