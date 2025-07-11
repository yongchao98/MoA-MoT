import math

def calculate_and_print_bound():
    """
    Calculates the constant for the upper bound on the systole
    and prints the final equation.
    """
    # Based on the theorem V >= 0.0883 * systole^3 by Colin Adams.
    # We rearrange this to get systole^3 <= (1/0.0883) * V.
    
    # The coefficient in the theorem.
    c_adams = 0.0883
    
    # The constant C in the inequality systole^3 <= C * V.
    C = 1 / c_adams
    
    # The original values in the rearranged inequality.
    numerator = 1
    denominator = c_adams

    # We will express the bound for k_{k,inf}^3.
    # The equation is: k_{k,inf}^3 <= C * V
    
    print("Based on a known theorem for hyperbolic 3-manifolds, the derived upper bound is:")
    print(f"k_k_inf**3 <= ({numerator}/{denominator}) * V")
    print("Which calculates to:")
    print(f"k_k_inf**3 <= {C:.4f} * V")


calculate_and_print_bound()