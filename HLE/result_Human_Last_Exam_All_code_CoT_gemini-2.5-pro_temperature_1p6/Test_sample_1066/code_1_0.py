import math

def print_lower_bound_formula():
    """
    This function prints the derived lower bound for the expected watermark score E[S].
    """
    
    # The lower bound for the expected score at a single position t is:
    # E[S_t] >= H(D_t) + (pi^2 / 6) - 1
    # Summing over all n tokens:
    # E[S] >= sum_{t=1 to n} (H(D_t) + (pi^2 / 6) - 1)
    # E[S] >= (sum_{t=1 to n} H(D_t)) + n * ((pi^2 / 6) - 1)
    # Using the definition of average entropy alpha = (1/n) * sum(H(D_t)),
    # we get sum(H(D_t)) = n * alpha.
    # Substituting this back gives the final lower bound formula.

    pi_squared_over_6 = math.pi**2 / 6
    one = 1.0

    print("The lower bound for E[S] is given by the formula:")
    # The prompt requires outputting each number in the final equation.
    # The formula is n * (alpha + pi^2/6 - 1)
    print(f"E[S] >= n * (alpha + {pi_squared_over_6:.4f} - {one:.4f})")
    print(f"E[S] >= n * (alpha + {pi_squared_over_6 - one:.4f})")
    print("\nIn symbolic form, the equation is:")
    print("E[S] >= n * (alpha + (pi^2 / 6) - 1)")

print_lower_bound_formula()
