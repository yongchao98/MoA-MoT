import math

def solve_immanant_problem(n: int):
    """
    This function solves the problem based on the step-by-step plan.
    
    The plan concludes that the matrix Mn can be constructed as a rank-1 nilpotent matrix
    to satisfy the problem's constraints under a reasonable interpretation. The specific
    matrix Mn is built from vectors r = (1, 1, ..., 1) and c = (1, 1, ..., 1, -(n-1)).

    For this matrix, it was shown that all immanants are zero except for the permanent.
    The permanent is calculated as -(n-1) * n!.
    
    Since n must be an integer >= 2 for an n x n non-trivial nilpotent matrix,
    the permanent -(n-1) * n! is always negative.
    The set of immanant values is {0, -(n-1) * n!}.
    The largest value in this set is 0.
    
    The final logic is independent of the value of n (for n>=2).
    """

    if n < 2:
        # A 1x1 nilpotent matrix with a non-zero entry is not possible.
        # The problem implies n >= 2.
        print("The problem is ill-defined for n < 2.")
        return

    # Based on the reasoning, the matrix M_n can be chosen such that
    # its permanent is negative. For a rank-1 matrix, all other
    # irreducible character immanants are zero.
    
    # The set of immanants is {0, Permanent}.
    # Permanent = -(n-1) * n!
    # For n >= 2, Permanent is negative.
    
    largest_immanant = 0
    
    # The problem asks to output the final equation.
    # The permanent is the only non-zero immanant, with value -(n-1)*n!
    # All other immanants are 0.
    # The largest immanant is max(0, -(n-1)*n!)
    
    perm_val = -(n - 1) * math.factorial(n)
    
    # We print the components of the comparison that leads to the answer.
    print(f"The non-zero immanant (permanent) is: {perm_val}")
    print(f"All other immanants are: 0")
    print(f"The largest immanant is max(0, {perm_val})")
    print(f"Final Answer: {largest_immanant}")


# The user did not specify a value for n, so we can demonstrate with a common choice like n=4.
# The logic holds for any n >= 2, and the final answer for the largest immanant is always 0.
n = 4 
solve_immanant_problem(n)