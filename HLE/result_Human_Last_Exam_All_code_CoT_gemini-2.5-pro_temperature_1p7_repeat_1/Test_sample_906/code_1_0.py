def solve_steady_state_probability():
    """
    This function explains the derivation of the steady-state probability pi_0 for
    a specific birth-death Markov process.
    The rates are:
    - Birth rate (i -> i+1): lambda_i = lambda / (i + 1)
    - Death rate (i -> i-1): mu_i = mu
    The solution is expressed in terms of rho = lambda / mu.
    """
    # Define symbols for printing
    rho = "ρ"
    pi_0 = "π₀"
    pi_n = "πₙ"
    
    print("Step-by-step derivation of the steady-state probability π₀:")
    print("-" * 55)
    
    # Step 1: Explain the derivation of π_n in terms of π_0
    print("1. For a birth-death process, the steady-state probabilities satisfy the detailed balance equations:")
    print("   πᵢ * λᵢ = πᵢ₊₁ * μᵢ₊₁ which gives πᵢ₊₁ = πᵢ * (λᵢ / μᵢ₊₁)")
    print(f"   Given λᵢ = λ/(i+1) and μᵢ = μ, we can find πₙ in terms of {pi_0}.")
    print(f"   Let {rho} = λ/μ.")
    print(f"   For n=1: π₁ = π₀ * (λ₀/μ₁) = π₀ * (λ/μ) = {pi_0} * {rho}¹/1!")
    print(f"   For n=2: π₂ = π₁ * (λ₁/μ₂) = ({pi_0}*{rho}) * ( (λ/2)/μ ) = {pi_0} * {rho}²/2!")
    print(f"   For n=3: π₃ = π₂ * (λ₂/μ₃) = ({pi_0}*{rho}²/2!) * ( (λ/3)/μ ) = {pi_0} * {rho}³/3!")
    print(f"   The general formula is: {pi_n} = {pi_0} * ({rho}ⁿ / n!)")
    print("-" * 55)

    # Step 2: Use the normalization condition
    print("2. The sum of all probabilities must be 1: Σπₙ = 1 (from n=0 to ∞).")
    print(f"   Substituting our formula: Σ [{pi_0} * ({rho}ⁿ / n!)] = 1")
    print(f"   Factoring out {pi_0}: {pi_0} * Σ({rho}ⁿ / n!) = 1")
    print("-" * 55)

    # Step 3: Identify the series and solve for π_0
    print("3. The series Σ({rho}ⁿ/n!) is the Taylor series for the exponential function e^{rho}.")
    print(f"   The equation becomes: {pi_0} * e^{rho} = 1.")
    print(f"   Solving for {pi_0}, we get: {pi_0} = 1 / e^{rho}")
    print("-" * 55)
    
    # Final Answer
    print("The final equation for the steady-state probability π₀ is:")
    # We output the '-1' to show all numbers in the final equation.
    print(f"π₀ = e ** (-1 * ρ)")

# Execute the derivation
solve_steady_state_probability()