import math

def calculate_minimum_curvature_cost(n, d):
    """
    Calculates and explains the minimum curvature cost for an NGD update.

    Args:
        n (int): The number of training samples.
        d (int): The dimension of the square weight matrix.
    """
    if not n < d:
        print(f"Warning: The condition n < d is not met (n={n}, d={d}). The analysis assumes n < d.")
        # We can still proceed to show the formulas
    
    # Number of parameters
    p = d * d

    print("Analysis of NGD Curvature Cost\n" + "="*35)
    print(f"Given parameters: n = {n} (samples), d = {d} (layer dimension)")
    print(f"The network has a single d x d layer, so the number of parameters is p = d*d = {p}.\n")

    # Step 1: Naive Approach
    print("1. Naive Inversion Cost")
    print("-" * 25)
    print(f"The NGD update requires inverting the matrix (F + alpha*I), where F is the Fisher matrix.")
    print(f"This matrix is of size p x p, which is {p} x {p}.")
    print(f"A direct inversion of a matrix of this size has a computational cost of O(p^3).")
    print(f"Substituting p = d^2, the naive cost is O((d^2)^3) = O(d^6).")
    try:
        # Use scientific notation for very large numbers
        naive_cost_exponent = d**6
        print(f"For d={d}, this complexity is O({d}^6) = O({'{:e}'.format(naive_cost_exponent)}), which is computationally prohibitive.\n")
    except OverflowError:
        print(f"For d={d}, this complexity is O({d}^6), which is computationally prohibitive.\n")


    # Step 2: Efficient Approach using Woodbury Identity
    print("2. Efficient Inversion Cost")
    print("-" * 25)
    print("The Fisher matrix F for a least-squares loss has a low-rank structure: F = J^T * J.")
    print(f"Here, J is the n x p (i.e., {n} x {p}) Jacobian matrix.")
    print("Since n < p, the rank of F is at most n.")
    print("\nWe can exploit this using the Woodbury Matrix Identity to avoid the large inversion.")
    print("This identity transforms the problem from inverting a p x p matrix to inverting an n x n matrix.")
    print(f"The smaller matrix to be inverted is of size n x n, which is {n} x {n}.")
    
    # Step 3: Final Minimum Cost
    print("\n3. Minimum Achievable Cost")
    print("-" * 25)
    print(f"The cost of inverting this smaller {n}x{n} matrix is O(n^3).")
    efficient_cost_val = n**3
    print(f"For n={n}, this complexity is O({n}^3) = O({efficient_cost_val}).")

    print("\n" + "="*35)
    print("The minimum achievable cost for the inversion operation in the NGD update is therefore O(n^3).")
    print("Final Equation: Cost = O(n^k)")
    print("In this equation, the base is n (the number of samples) and the exponent k is 3.")

# Example values as per the problem description (n < d)
n_samples = 50
d_dimension = 200

calculate_minimum_curvature_cost(n_samples, d_dimension)