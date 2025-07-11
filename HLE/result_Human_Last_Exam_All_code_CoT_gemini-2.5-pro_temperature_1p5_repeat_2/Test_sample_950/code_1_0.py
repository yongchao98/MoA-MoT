import math

def solve_grassmannian_torsion_rank():
    """
    Calculates the 2-rank of the integral cohomology of the real Grassmannian Gr_k(R^n).
    This corresponds to the rank of the torsion subgroup.
    """
    # Parameters for the Grassmannian Gr_k(R^n) from the problem
    # space of 3-subspaces of R^5 means k=3, n=5.
    k = 3
    n = 5

    print(f"Calculating the 2-rank of the torsion subgroup for the integral cohomology of Gr_({k})(R^{n}).")
    print(f"This is the number of integers i (from 1 to {k}) for which the binomial coefficient C(n-k-1+i, i) is odd.")
    print("-" * 30)

    rank_count = 0
    summands = []

    for i in range(1, k + 1):
        # Formula for the binomial coefficient's parameters
        N = n - k - 1 + i
        K = i
        
        # Calculate the binomial coefficient C(N, K)
        try:
            # math.comb requires Python 3.8+
            c_val = math.comb(N, K)
            # Check if the result is odd
            is_odd = (c_val % 2 != 0)
        except ValueError:
            # This can happen if K > N, in which case C(N,K)=0
            c_val = 0
            is_odd = False
            
        term_contribution = 1 if is_odd else 0
        summands.append(str(term_contribution))
        rank_count += term_contribution
        
        parity_string = "odd" if is_odd else "even"
        
        print(f"For i = {i}:")
        print(f"  N = {n}-{k}-1+{i} = {N}")
        print(f"  K = {i}")
        print(f"  C({N}, {K}) = {c_val}, which is {parity_string}.")
        print(f"  This contributes {term_contribution} to the total rank.")
        print()

    # Create the final equation string like "0 + 1 + 0"
    equation_string = " + ".join(summands)

    print("-" * 30)
    print("The total rank is the sum of these contributions:")
    print(f"Final Equation: {equation_string} = {rank_count}")

solve_grassmannian_torsion_rank()