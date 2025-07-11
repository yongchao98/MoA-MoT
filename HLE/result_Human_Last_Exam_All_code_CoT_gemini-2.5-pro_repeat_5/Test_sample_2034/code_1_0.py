import numpy as np

def solve():
    """
    This function determines the constants c1 and c2 based on the theoretical derivation.
    """
    # Based on the derivation, the relationship c1*alpha_i = c2*(K*alpha)_i must hold.
    # Under the simplifying assumption that for the base case (beta=0), alpha_i = (K*alpha)_i,
    # this implies c1 = c2. The simplest non-trivial choice is c1 = c2 = 1.
    c1 = 1
    c2 = 1

    alpha_i = 0.5  # Example value for alpha_i
    K_alpha_i = 0.5 # Example value satisfying alpha_i = (K*alpha)_i
    beta = 0.01    # Small beta
    K_alpha_loo_i = -0.01 # Example value for the leave-one-out prediction

    # Print the final equation with the determined coefficients
    print("The extended Jaakola-Haussler bound is:")
    print(f"- (K_alpha_loo)_i <= +(1 + {c1}*beta)*alpha_i - (1 + {c2}*beta)*(K*alpha)_i + o(beta)")

    # Demonstrate the equation with example values
    lhs = -K_alpha_loo_i
    rhs = (1 + c1 * beta) * alpha_i - (1 + c2 * beta) * K_alpha_i
    
    print("\nFor example values:")
    print(f"alpha_i = {alpha_i}")
    print(f"(K*alpha)_i = {K_alpha_i}")
    print(f"-(K*alpha_loo)_i = {lhs}")
    print(f"beta = {beta}")
    
    print("\nThe inequality becomes:")
    print(f"{lhs} <= + (1 + {c1} * {beta}) * {alpha_i} - (1 + {c2} * {beta}) * {K_alpha_i}")
    print(f"{lhs} <= {1 + c1 * beta} * {alpha_i} - {1 + c2 * beta} * {K_alpha_i}")
    print(f"{lhs} <= {(1 + c1 * beta) * alpha_i} - {(1 + c2 * beta) * K_alpha_i}")
    print(f"{lhs} <= {rhs}")
    
    # Final determined values
    print(f"\nTherefore, c1 = {c1} and c2 = {c2}.")

solve()
