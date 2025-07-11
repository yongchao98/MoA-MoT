import math

def solve_dynamics_problem():
    """
    Calculates the number of iterations to distinguish a chaotic orbit
    and performs the final required computation.
    """

    # Step 1: Define the given parameters from the problem.
    # n represents the period, indicating a chaotic regime.
    n = 3
    # p is the precision in significant digits for Era B.
    p = 7
    # The divisor specified in the final calculation step.
    divisor = 12

    # Step 2: Define the constants for the calculation.
    # The Lyapunov exponent (λ) for the logistic map at r=4 (canonical chaos) is ln(2).
    lambda_exponent = math.log(2)
    # The natural logarithm of 10 is needed for the conversion from base-10 precision.
    ln_10 = math.log(10)

    # Step 3: Calculate T(n,p), the predictability horizon.
    # The formula is T = p * ln(10) / λ.
    # This can also be written as T = p * log_2(10).
    T_np = p * ln_10 / lambda_exponent

    # Step 4: Calculate the final result required by the problem.
    final_result = math.ceil(T_np / divisor)

    # Step 5: Print the explanation and the final equation with all its numbers.
    print("--- Calculation Steps ---")
    print("Problem: Calculate ceil(T(n,p) / 12)")
    print(f"Given parameters: n = {n}, p = {p}")
    print("\nFormula for predictability horizon T(n,p):")
    print("T(n,p) = p * ln(10) / λ")
    print(f"Using λ = ln(2) for canonical chaos (logistic map with r=4).")
    
    print("\nSubstituting the values into the formula:")
    print(f"T({n},{p}) = {p} * ln(10) / ln(2)")
    print(f"T({n},{p}) = {p} * {ln_10:.5f} / {lambda_exponent:.5f}")
    print(f"T({n},{p}) ≈ {T_np:.4f}")

    print("\nFinal computation:")
    print(f"ceil(T({n},{p}) / {divisor}) = ceil({T_np:.4f} / {divisor})")
    print(f"Result = {final_result}")

solve_dynamics_problem()