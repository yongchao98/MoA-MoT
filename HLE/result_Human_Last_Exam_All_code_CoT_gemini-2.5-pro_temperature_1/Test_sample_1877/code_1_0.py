import math

def solve_problem():
    """
    Calculates T(n,p) and the final result based on the shadowing time of unstable periodic orbits.
    """
    # Parameters from the problem
    n = 3  # Period of the orbit
    p = 7  # Era B precision (7 significant digits)

    # For the logistic map in a chaotic regime (e.g., r=4), the stability
    # multiplier for the unstable period-3 orbits is a known value.
    Lambda_n = 8.0

    # Calculate ln(10) and ln(Lambda_n)
    ln_10 = math.log(10)
    ln_Lambda_n = math.log(Lambda_n)

    # Calculate T(n,p) using the formula: T = (n * p * ln(10)) / ln(Lambda_n)
    T_np = (n * p * ln_10) / ln_Lambda_n

    # The final step is to calculate ceil(T(n,p) / 12)
    final_value = 12.0
    result = math.ceil(T_np / final_value)

    # Output the steps of the calculation
    print("The formula for the shadowing time T(n,p) is:")
    print("T(n, p) = (n * p * ln(10)) / ln(Λ_n)\n")
    print(f"Given parameters:")
    print(f"n = {n}")
    print(f"p = {p}")
    print(f"Λ_n = {Lambda_n} (Stability multiplier for unstable period-3 orbit at r=4)\n")
    print(f"Intermediate calculation for T({n},{p}):")
    print(f"T({n},{p}) = ({n} * {p} * {ln_10:.6f}) / {ln_Lambda_n:.6f} = {T_np:.6f}\n")
    print("Final calculation as required by the problem:")
    print(f"ceil(T({n},{p}) / {final_value}) = ceil({T_np:.6f} / {final_value}) = ceil({T_np/final_value:.6f}) = {int(result)}")
    print("\n<<<{}>>>".format(int(result)))


solve_problem()