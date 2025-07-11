from fractions import Fraction

def solve_variance():
    """
    Calculates the variance of Y based on the order statistics method.
    """
    n = 4

    # Step 1: Determine the probability that Y equals a specific order statistic X_(k).
    # This is derived from case analysis:
    # - If X1 = X_(1), Y = X_(3). (Prob 1/4)
    # - If X1 = X_(4), Y = X_(2). (Prob 1/4)
    # - If X1 = X_(2), Y is X_(1) or X_(3) each with prob 1/2. (Total prob 1/4)
    # - If X1 = X_(3), Y is X_(2) or X_(4) each with prob 1/2. (Total prob 1/4)
    
    # P(Y = X_(k))
    P_Y_is_Xk = {
        1: Fraction(1, 4) * Fraction(1, 2),  # Case X1=X_(2) and d1 > d2
        2: Fraction(1, 4) + Fraction(1, 4) * Fraction(1, 2), # Case X1=X_(4) + Case X1=X_(3) and d2 > d3
        3: Fraction(1, 4) + Fraction(1, 4) * Fraction(1, 2), # Case X1=X_(1) + Case X1=X_(2) and d1 < d2
        4: Fraction(1, 4) * Fraction(1, 2),  # Case X1=X_(3) and d2 < d3
    }

    print("Probabilities P(Y = X_(k)):")
    for k in range(1, n + 1):
        print(f"  P(Y = X_({k})) = {P_Y_is_Xk[k]}")

    # Step 2: Calculate the first and second moments of the order statistics X_(k) for n=4.
    # E[X_(k)] = k / (n+1)
    # E[X_(k)^2] = k*(k+1) / ((n+1)*(n+2))
    E_Xk = {k: Fraction(k, n + 1) for k in range(1, n + 1)}
    E_Xk_sq = {k: Fraction(k * (k + 1), (n + 1) * (n + 2)) for k in range(1, n + 1)}
    
    print("\nMoments of order statistics X_(k):")
    for k in range(1, n + 1):
        print(f"  E[X_({k})] = {E_Xk[k]}, E[X_({k})^2] = {E_Xk_sq[k]}")

    # Step 3: Calculate E[Y] and E[Y^2] using the law of total expectation.
    E_Y = sum(P_Y_is_Xk[k] * E_Xk[k] for k in range(1, n + 1))
    E_Y_sq = sum(P_Y_is_Xk[k] * E_Xk_sq[k] for k in range(1, n + 1))
    
    print(f"\nExpected value E[Y] = {E_Y}")
    print(f"Expected value of square E[Y^2] = {E_Y_sq}")

    # Step 4: Calculate the variance of Y.
    # Var(Y) = E[Y^2] - (E[Y])^2
    var_Y = E_Y_sq - E_Y**2
    
    print("\nCalculating the variance:")
    print(f"Var(Y) = E[Y^2] - (E[Y])^2")
    print(f"Var(Y) = {E_Y_sq} - ({E_Y})^2")
    print(f"Var(Y) = {E_Y_sq} - {E_Y**2}")
    print(f"Var(Y) = {var_Y}")
    
solve_variance()