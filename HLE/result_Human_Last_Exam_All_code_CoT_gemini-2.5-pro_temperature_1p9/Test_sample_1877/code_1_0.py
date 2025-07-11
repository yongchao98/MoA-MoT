import math

def solve_problem():
    """
    Calculates T(n,p) for n=3 and p=7 and derives the final answer.
    """
    # Step 1: Define the parameters based on the problem statement
    # n = 3 (period-3 orbit comparison)
    # p = 7 (Era B: 7 significant digits)
    p = 7
    n = 3

    # Step 2: Determine the Lyapunov exponent (λ).
    # To find the *minimum* number of iterations T, we need the *maximum* rate of
    # divergence, which corresponds to the maximum Lyapunov exponent for the logistic map.
    # This occurs at r=4, where λ = ln(2).
    lambda_exp = math.log(2)

    # Step 3: Use the formula T(n,p) = (p * ln(10)) / λ
    # This formula calculates the number of iterations for an initial error of 10^(-p)
    # to grow to order 1.
    ln_10 = math.log(10)
    T_3_7 = (p * ln_10) / lambda_exp

    # Step 4: Calculate the final value as requested by the problem.
    # Final value = ceil(T(3,7) / 12)
    denominator = 12
    final_value_float = T_3_7 / denominator
    final_answer = math.ceil(final_value_float)

    # Step 5: Print the explanation and each number in the final equations.
    print("The formula for the predictability horizon T(n, p) is:")
    print("T(n, p) = (p * ln(10)) / λ\n")

    print("For this problem, the parameters are:")
    print(f"p = {p} (from 'Era B' precision)")
    print(f"n = {n}")
    print(f"λ = ln(2) (maximal Lyapunov exponent for fastest divergence)")
    print("-" * 25)

    print("Calculation for T(3, 7):")
    print(f"T(3, 7) = ({p} * ln(10)) / ln(2)")
    print(f"T(3, 7) = ({p} * {ln_10:.5f}) / {lambda_exp:.5f}")
    print(f"T(3, 7) = {p * ln_10:.5f} / {lambda_exp:.5f}")
    print(f"T(3, 7) ≈ {T_3_7:.5f}\n")

    print("Final calculation:")
    print(f"ceil(T(3, 7) / {denominator}) = ceil({T_3_7:.5f} / {denominator})")
    print(f"                      = ceil({final_value_float:.5f})")
    print(f"                      = {final_answer}")


solve_problem()