import math

def print_lower_bound_formula():
    """
    This function prints the derived lower bound for the expected detection score E[S].
    """
    n_str = "n"
    alpha_str = "alpha"
    gamma_str = "gamma"
    
    # The derived lower bound for E[S] is n * (alpha + gamma).
    # We will print this formula, including the numerical value of gamma.
    
    gamma_val = math.gamma(1) # This is not gamma, this is Gamma(1)=0! = 1
    # The euler_mascheroni constant is not in the standard math library.
    # We define it manually.
    euler_mascheroni_gamma = 0.5772156649015328
    
    # The problem statement requires printing the numbers in the final equation.
    # The only numerical constant in our derived bound is gamma.
    
    print("The derived lower bound for E[S] is:")
    print(f"E[S] >= {n_str} * ({alpha_str} + {gamma_str})")
    print("\nSubstituting the approximate value for the Euler-Mascheroni constant (gamma):")
    print(f"E[S] >= {n_str} * ({alpha_str} + {euler_mascheroni_gamma})")

print_lower_bound_formula()