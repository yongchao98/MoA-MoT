import numpy as np

def solve_markov_chain():
    """
    Solves for the steady-state probabilities of the given Kolmogorov-Chepmen system
    and calculates P0(inf) + P1(inf).
    """
    # Given lambda values
    lambda_01 = 0.019
    lambda_10 = 0.65
    lambda_12 = 0.4
    lambda_21 = 0.392
    lambda_23 = 0.008
    lambda_31 = 0.008

    print("Step 1: Set up the steady-state equations (P'(t) = 0). Let p_i = P_i(inf).")
    print(f"Eq 1: -{lambda_01}*p0 + {lambda_10}*p1 = 0")
    print(f"Eq 2: {lambda_12}*p1 - ({lambda_21} + {lambda_23})*p2 = 0")
    print(f"Eq 3: {lambda_23}*p2 - {lambda_31}*p3 = 0")
    print("Eq 4: p0 + p1 + p2 + p3 = 1\n")

    print("Step 2: Solve for p0, p2, p3 in terms of p1.")
    # From Eq 3: lambda_23 * p2 = lambda_31 * p3
    # 0.008 * p2 = 0.008 * p3  => p3 = p2
    print(f"From Eq 3: {lambda_23}*p2 = {lambda_31}*p3  => p3 = ({lambda_23}/{lambda_31})*p2 = {lambda_23/lambda_31}*p2. So, p3 = p2.")

    # From Eq 2: lambda_12 * p1 = (lambda_21 + lambda_23) * p2
    lambda_21_23_sum = lambda_21 + lambda_23
    # 0.4 * p1 = (0.392 + 0.008) * p2 = 0.4 * p2 => p2 = p1
    print(f"From Eq 2: {lambda_12}*p1 = ({lambda_21} + {lambda_23})*p2 = {lambda_21_23_sum}*p2  => p2 = ({lambda_12}/{lambda_21_23_sum})*p1 = {lambda_12/lambda_21_23_sum}*p1. So, p2 = p1.")
    print("Therefore, p1 = p2 = p3.\n")

    # From Eq 1: lambda_10 * p1 = lambda_01 * p0
    c0 = lambda_10 / lambda_01
    print(f"From Eq 1: {lambda_10}*p1 = {lambda_01}*p0 => p0 = ({lambda_10}/{lambda_01})*p1 = {c0:.4f}*p1.\n")
    
    print("Step 3: Substitute these into the normalization equation p0 + p1 + p2 + p3 = 1.")
    # p0 + p1 + p2 + p3 = 1
    # (c0 * p1) + p1 + p1 + p1 = 1
    # (c0 + 3) * p1 = 1
    total_coeff_p1 = c0 + 3
    print(f"({c0:.4f})*p1 + p1 + p1 + p1 = 1")
    print(f"({c0:.4f} + 3)*p1 = 1")
    print(f"{total_coeff_p1:.4f}*p1 = 1\n")

    print("Step 4: Solve for p1 and p0.")
    p1 = 1 / total_coeff_p1
    p0 = c0 * p1
    print(f"p1 = 1 / {total_coeff_p1:.4f} = {p1:.6f}")
    print(f"p0 = {c0:.4f} * p1 = {c0:.4f} * {p1:.6f} = {p0:.6f}\n")

    print("Step 5: Calculate the final sum P0(inf) + P1(inf) = p0 + p1.")
    result = p0 + p1
    print(f"The final sum is p0 + p1 = {p0:.6f} + {p1:.6f} = {result:.6f}")
    
    # Alternative calculation for display
    # p0+p1 = (c0+1)*p1 = (c0+1)/(c0+3)
    # (0.65+0.019)/0.707 = 0.669 / 0.707
    print("\n--- Final Equation ---")
    print(f"P0(+inf) + P1(+inf) = p0 + p1")
    print(f"p0 + p1 = (({lambda_10}/{lambda_01}) * p1) + p1 = (({lambda_10}/{lambda_01}) + 1) * p1")
    print(f"p1 = 1 / (({lambda_10}/{lambda_01}) + 1 + 1 + 1) = 1 / (({lambda_10}/{lambda_01}) + 3)")
    print(f"p0 + p1 = (({lambda_10}/{lambda_01}) + 1) / (({lambda_10}/{lambda_01}) + 3)")
    print(f"p0 + p1 = (({lambda_10} + {lambda_01})/{lambda_01}) / (({lambda_10} + 3 * {lambda_01})/{lambda_01})")
    print(f"p0 + p1 = ({lambda_10} + {lambda_01}) / ({lambda_10} + 3 * {lambda_01})")
    numerator = lambda_10 + lambda_01
    denominator = lambda_10 + 3 * lambda_01
    print(f"p0 + p1 = {numerator} / {denominator}")
    final_value = numerator / denominator
    print(f"p0 + p1 = {final_value}")

solve_markov_chain()