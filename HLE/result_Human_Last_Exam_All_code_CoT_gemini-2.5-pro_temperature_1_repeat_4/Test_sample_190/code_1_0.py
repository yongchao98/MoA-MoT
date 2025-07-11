import sympy as sp

def solve_markov_transience():
    """
    This function symbolically calculates the critical value of 'c' for which
    the given Markov chain becomes transient.
    """
    # Define symbols for state k and parameter c
    k = sp.Symbol('k', real=True, positive=True)
    c = sp.Symbol('c', real=True)

    # Define the transition probabilities for a large state k
    P_km2 = sp.Rational(1, 4)
    P_kp2 = sp.Rational(1, 4)
    P_km1 = sp.Rational(1, 4) - c/k
    P_kp1 = sp.Rational(1, 4) + c/k

    # We use a Lyapunov function f(x) = log(x) to test for transience.
    # We need to compute the expected drift: E[log(X_{n+1}) - log(X_n) | X_n = k]
    # This is equivalent to E[log(X_{n+1}/k) | X_n = k]
    expected_log_ratio = (
        P_km2 * sp.log(1 - 2/k) +
        P_kp2 * sp.log(1 + 2/k) +
        P_km1 * sp.log(1 - 1/k) +
        P_kp1 * sp.log(1 + 1/k)
    )

    # To analyze the behavior for large k, we perform a Taylor series expansion
    # with respect to 1/k around 0. We need to go to the (1/k)^2 term
    # as the (1/k) terms cancel out.
    series_expansion = sp.series(expected_log_ratio, k, sp.oo, 3).removeO()

    # The sign of the drift for large k is determined by the coefficient
    # of the dominant term in the expansion, which is the 1/k^2 term.
    # We isolate this coefficient.
    drift_coefficient = sp.simplify(series_expansion * k**2)

    # For the chain to be transient, the drift must be positive.
    # This means the drift_coefficient must be positive.
    # The critical value is where this coefficient is zero.
    critical_equation = sp.Eq(drift_coefficient, 0)
    
    # Solve the equation for c to find the critical value.
    solution = sp.solve(critical_equation, c)
    infimum_c = solution[0]

    # Print the step-by-step derivation
    print("Plan: Analyze the drift of f(X_n) = log(X_n) to find the condition for transience.")
    print("\n1. The expected one-step change for log(X_n) from state k is:")
    print(f"   Drift = E[log(X_{n+1}/k) | X_n = k]")
    
    print("\n2. For large k, the Taylor series expansion of this drift is:")
    print(f"   Drift ≈ {series_expansion}")
    
    print("\n3. For the chain to be transient, this drift must be positive.")
    print("   This requires the coefficient of the dominant term to be positive.")
    print(f"   The coefficient is: {drift_coefficient}")
    
    print("\n4. The condition for transience is: " + str(drift_coefficient > 0))
    
    print("\n5. The boundary between recurrence and transience is where the coefficient is zero.")
    print("   The final equation for this critical point is:")
    c_coeff = drift_coefficient.coeff(c)
    const_term = drift_coefficient.subs(c, 0)
    print(f"   {int(c_coeff)}*c + ({const_term}) = 0")
    
    print("\n6. Solving for c, the infimum of the set of c for which the chain is transient is:")
    print(f"   c = {infimum_c}")

solve_markov_transience()