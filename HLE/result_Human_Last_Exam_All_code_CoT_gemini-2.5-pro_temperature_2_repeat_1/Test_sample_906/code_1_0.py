def solve_steady_state_probability():
    """
    This function derives the steady-state probability pi_0 for the given
    birth-death process and prints the derivation step-by-step.
    """
    # Define symbolic representations for the variables
    lambda_rate = "λ"
    mu_rate = "μ"
    rho = "ρ"
    pi_0 = "π_0"
    pi_n = "π_n"

    print("1. The general formula for π_n in terms of π_0 for a birth-death process is:")
    print(f"   {pi_n} = {pi_0} * Π_{{i=0 to n-1}} (λ_i / μ_{{i+1}})")
    print()

    print("2. Substitute the given rates:")
    print(f"   Birth rate λ_i = {lambda_rate} / (i + 1)")
    print(f"   Death rate μ_i = {mu_rate}")
    print("   So, μ_{i+1} is always μ.")
    print()

    print("3. Derive the expression for π_n:")
    print(f"   {pi_n} = {pi_0} * Π_{{i=0 to n-1}} (({lambda_rate} / (i + 1)) / {mu_rate})")
    print(f"   {pi_n} = {pi_0} * [({lambda_rate}/(1*{mu_rate})) * ({lambda_rate}/(2*{mu_rate})) * ... * ({lambda_rate}/(n*{mu_rate}))]")
    print(f"   {pi_n} = {pi_0} * ({lambda_rate}^n / ({mu_rate}^n * n!))")
    print()

    print(f"4. Substitute ρ = {lambda_rate} / {mu_rate}:")
    print(f"   {pi_n} = {pi_0} * ({rho}^n / n!)")
    print()

    print("5. Use the normalization condition Σ π_n = 1:")
    print(f"   Σ_{{n=0 to ∞}} {pi_0} * ({rho}^n / n!) = 1")
    print(f"   {pi_0} * Σ_{{n=0 to ∞}} ({rho}^n / n!) = 1")
    print()
    
    print("6. Recognize the sum as the Taylor series for e^ρ:")
    print("   Σ_{n=0 to ∞} (ρ^n / n!) = e^ρ")
    print()

    print("7. Solve for π_0:")
    print(f"   {pi_0} * e^{rho} = 1")
    print(f"   {pi_0} = 1 / e^{rho}")
    print(f"   {pi_0} = e^(-{rho})")
    print()
    
    final_answer = "e^(-ρ)"
    print(f"The final expression for π_0 is: {final_answer}")

solve_steady_state_probability()