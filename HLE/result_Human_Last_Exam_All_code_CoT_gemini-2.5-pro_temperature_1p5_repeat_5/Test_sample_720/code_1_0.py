def get_minimum_curvature_cost():
    """
    This function analyzes the computational cost of the NGD update and
    prints the formula for the minimum achievable curvature cost.
    
    The cost is derived for a single-layer fully-connected network of size d x d,
    trained on n samples where n < d.
    
    The final cost is expressed in terms of d and n.
    """
    
    # The computational cost is determined by the most efficient algorithm
    # to compute the update. The algorithm involves a series of matrix operations.
    # Let's denote the cost of each step using d and n.
    
    # The cost is a polynomial in d and n. The coefficients represent the
    # approximate number of times a certain operation is performed.
    
    # Cost formula: c1 * d^2 * n + c2 * d * n^2 + c3 * n^3
    # From our step-by-step analysis, the coefficients are:
    c1 = 2
    c2 = 2
    c3 = 1

    # The detailed cost formula represents the sum of costs of the main
    # computational steps:
    # 2*(d^2*n) for matrix products involving d x d and d x n matrices
    # 2*(d*n^2) for matrix products involving d x n and n x n matrices and Gram matrix calculation
    # 1*(n^3)   for inverting an n x n matrix
    
    print("The minimum curvature cost is the cost of the most efficient algorithm to compute the NGD update.")
    print("The detailed cost formula, based on the number of floating point operations, is approximately:")
    print(f"Cost = {c1}*d^2*n + {c2}*d*n^2 + {c3}*n^3")
    
    # Given the condition n < d, the term d^2*n is asymptotically dominant.
    dominant_term = "O(d^2 * n)"
    
    print(f"\nSince n < d, the asymptotic complexity (Big-O notation) is determined by the largest term, which is:")
    print(dominant_term)


if __name__ == '__main__':
    get_minimum_curvature_cost()
