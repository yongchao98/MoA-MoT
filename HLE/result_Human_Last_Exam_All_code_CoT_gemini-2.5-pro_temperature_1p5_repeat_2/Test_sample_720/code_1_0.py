import sys

def calculate_minimum_curvature_cost(d, n):
    """
    Calculates and explains the minimum curvature cost for the given NGD update scenario.

    The cost is the total computational complexity of the inversion operation
    (F + alpha*I)^-1 * g.

    Args:
        d (int): The dimension of the neural network's single square layer.
        n (int): The number of training samples, where n < d.
    """
    if n >= d:
        print("Error: The analysis assumes n < d, but received d={} and n={}".format(d, n), file=sys.stderr)
        return

    num_params = d**2
    print("--- NGD Curvature Cost Analysis ---")
    print(f"Network Layer Size: {d} x {d}")
    print(f"Number of Parameters (p): d^2 = {num_params}")
    print(f"Number of Samples (n): {n}")
    print(f"Constraint: n < d ({n} < {d}) is met.\n")

    print("1. Naive Inversion Cost:")
    print(f"The Fisher matrix F is {num_params}x{num_params}. A direct inversion would cost O(p^3) = O((d^2)^3) = O(d^6).\n")

    print("2. Exploiting Low-Rank Structure:")
    print(f"The empirical Fisher matrix F is computed from n samples, so its rank is at most n.")
    print(f"Since n ({n}) is much smaller than the matrix dimension d^2 ({num_params}), F is a low-rank matrix.")
    print("We can use the Woodbury matrix identity to convert the large inversion into a smaller n x n inversion.\n")

    print("3. Cost of the Efficient Method:")
    cost_formation_str = "O(n^2 * d^2)"
    cost_inversion_str = "O(n^3)"
    cost_matvec_str = "O(n * d^2)"
    print(f"The cost is composed of three main parts:")
    print(f"  a) Forming an n x n matrix (K = JJ^T):  {cost_formation_str}")
    print(f"  b) Inverting the n x n matrix:          {cost_inversion_str}")
    print(f"  c) Matrix-vector products:             {cost_matvec_str}\n")

    print("4. Minimum Curvature Cost:")
    print(f"The total cost is the sum of these parts. Given n < d, the term {cost_formation_str} is the largest.")
    print("Therefore, the minimum achievable curvature cost is dominated by the cost of forming the n x n matrix.")
    
    final_cost_formula = "n^2 * d^2"
    print(f"\nFinal Dominant Cost Formula: {final_cost_formula}")
    
    # As requested, outputting each number in the final equation using example values
    n_squared = n**2
    d_squared = d**2
    total_ops = n_squared * d_squared
    
    print("\nIllustrative Calculation for the Dominant Term:")
    print(f"Using n={n} and d={d}:")
    print(f"Cost is proportional to -> {n}^2 * {d}^2 = {n_squared} * {d_squared} = {total_ops}")


if __name__ == '__main__':
    # Example values satisfying n < d
    d_example = 100
    n_example = 50
    calculate_minimum_curvature_cost(d_example, n_example)