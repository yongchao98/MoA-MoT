import math

def calculate_minimum_curvature_cost(n, d):
    """
    Calculates the minimum achievable curvature cost for the NGD update.

    The cost is measured in approximate floating-point operations (FLOPs).
    A FLOP is considered one multiplication and one addition.

    Args:
        n (int): The number of samples.
        d (int): The dimension of the weight matrix (d x d).

    Returns:
        None. Prints the cost analysis.
    """
    if not n < d:
        print("This analysis assumes the number of samples n is less than the dimension d.")
        return

    print(f"Calculating minimum curvature cost for n = {n} and d = {d}:\n")
    print("The minimum cost is achieved by using the Woodbury identity to avoid a direct O(d^3) inversion.")
    print("The total cost is the sum of FLOPs for the following steps:\n")

    # Step 1: Cost of Z1 = X @ g_W  (n x d) @ (d x d)
    # Approx FLOPs = n * d * d
    cost_step1 = n * d**2
    print(f"1. Cost of matrix multiplication (X @ g_W):        {n} * {d}^2 = {cost_step1:,.0f} FLOPs")

    # Step 2: Cost of G = X @ X.T  (n x d) @ (d x n)
    # Approx FLOPs = n * d * n
    cost_step2 = n**2 * d
    print(f"2. Cost of matrix multiplication (X @ X.T):       {n}^2 * {d} = {cost_step2:,.0f} FLOPs")

    # Step 3: Cost of solving the n x n linear system for d vectors
    # Approx FLOPs = n^3 (for LU decomp) + d * n^2 (for d solves)
    cost_step3_lu = n**3
    cost_step3_solve = d * n**2
    cost_step3 = cost_step3_lu + cost_step3_solve
    print(f"3. Cost of solving the n x n linear system:     ({n}^3 + {d} * {n}^2) = {cost_step3:,.0f} FLOPs")

    # Step 4: Cost of Z3 = X.T @ Z2  (d x n) @ (n x d)
    # Approx FLOPs = d * n * d
    cost_step4 = d**2 * n
    print(f"4. Cost of matrix multiplication (X.T @ Z2):      {d}^2 * {n} = {cost_step4:,.0f} FLOPs")

    # Summing up the costs
    total_cost = cost_step1 + cost_step2 + cost_step3 + cost_step4
    print("\n---------------------------------------------------------")
    print(f"Total Minimum Curvature Cost (approximate) = {total_cost:,.0f} FLOPs")
    print("---------------------------------------------------------")

    # Determine dominant term for symbolic complexity
    # Since n < d, the term n*d^2 (which is the same as d^2*n) is the largest.
    print(f"\nFor n < d, the computational complexity is dominated by the O(nd^2) terms.")
    print("Therefore, the minimum achievable curvature cost has a complexity of O(n*d^2).")


if __name__ == '__main__':
    # Example values as specified in the problem, where n < d
    n_samples = 50
    dimension = 500
    calculate_minimum_curvature_cost(n_samples, dimension)
