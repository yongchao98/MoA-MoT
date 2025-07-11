def solve_entropy_problem():
    """
    This function calculates the maximal entropy based on the derived optimal construction.
    In the optimal case, the total entropy is the sum of the entropies of two
    independent variables, s1 and s2, whose entropies are maximized at 1.
    """
    # H(s1) is the entropy of s1. Its maximum value under the constraint H(s1) <= 1 is 1.
    H_s1 = 1
    
    # H(s2) is the entropy of s2. Its maximum value under the constraint H(s2) <= 1 is 1.
    H_s2 = 1
    
    # The maximal entropy H(x,y,z,s1,s2) is H(s1) + H(s2) in the optimal case.
    max_entropy = H_s1 + H_s2
    
    # Print the final equation with the numbers used.
    print(f"The maximal entropy is achieved when H(s1) and H(s2) are maximized.")
    print(f"The calculation is based on the identity H(x,y,z,s1,s2) = H(s1) + H(s2) for the optimal distribution.")
    print(f"Final Equation: {H_s1} + {H_s2} = {max_entropy}")

solve_entropy_problem()