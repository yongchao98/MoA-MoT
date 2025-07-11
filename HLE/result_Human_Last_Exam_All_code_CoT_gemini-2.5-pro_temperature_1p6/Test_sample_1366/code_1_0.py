import math

def solve_voa_problem():
    """
    Solves the VOA problem and prints the answer in the required format.
    """

    # Part (a): Theoretical result on VOA decomposition.
    # The decomposition is known to exist and is unique for irreducible modules.
    answer_a = "Yes; No"

    # Part (b): Top-level dimension from the problem's definition.
    # The top level of L(p)_n is rho_n, which is (n+1)-dimensional.
    answer_b = "n+1"

    # Part (c): Calculation of the minimal conformal weight for p=2.
    # The formula for the conformal weight is h_n = p*n*(n+2)/8.
    # We need the minimal non-zero weight, which occurs at n=1.
    p = 2
    n = 1

    # Numerator and denominator of the h_n formula
    numerator = p * n * (n + 2)
    denominator = 8

    # Simplify the resulting fraction
    common_divisor = math.gcd(numerator, denominator)
    simple_num = numerator // common_divisor
    simple_den = denominator // common_divisor
    
    # Format the answer for (c) to show the equation as requested.
    equation_str = f"({p}*{n}*({n}+2))/{denominator}"
    result_str = f"{simple_num}/{simple_den}"
    answer_c = f"{equation_str} = {result_str}"

    # Print the final combined answer.
    print(f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}")

solve_voa_problem()