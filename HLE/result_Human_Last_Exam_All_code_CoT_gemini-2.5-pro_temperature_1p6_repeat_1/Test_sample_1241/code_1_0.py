import sys
# It's better to use Fraction for precision
# But floating point should be sufficient here
# from fractions import Fraction

def solve_stationary_distribution():
    """
    Solves the stationary distribution for the given Kolmogorov-Chepmen system
    and calculates P0(inf) + P1(inf).
    """
    
    # --- Given parameters ---
    lambda_01 = 0.019
    lambda_10 = 0.65
    lambda_12 = 0.4
    lambda_21 = 0.392
    lambda_23 = 0.008
    lambda_31 = 0.008

    print("Step 1: Set up the steady-state equations (P_i' = 0).")
    print("Let P_i denote P_i(+inf).")
    print(f" (1) 0 = -{lambda_01}*P_0 + {lambda_10}*P_1")
    # We ignore the second equation as it contains a likely typo and is redundant.
    # The correct form would be P_1' = +lambda_01*P_0 - (lambda_10 + lambda_12)*P_1 + ...
    print(f" (2) 0 = {lambda_12}*P_1 - ({lambda_21} + {lambda_23})*P_2")
    print(f" (3) 0 = {lambda_23}*P_2 - {lambda_31}*P_3")
    print(" (4) P_0 + P_1 + P_2 + P_3 = 1")
    print("-" * 30)

    print("Step 2: Express P_0, P_2, and P_3 in terms of P_1.")
    # From (1): lambda_01 * P_0 = lambda_10 * P_1
    c0_ratio = lambda_10 / lambda_01
    print(f"From (1): P_0 = ({lambda_10} / {lambda_01}) * P_1 = {c0_ratio:.4f} * P_1")

    # From (2): (lambda_21 + lambda_23) * P_2 = lambda_12 * P_1
    l21_plus_l23 = lambda_21 + lambda_23
    c2_ratio = lambda_12 / l21_plus_l23
    print(f"From (2): P_2 = ({lambda_12} / ({lambda_21} + {lambda_23})) * P_1 = ({lambda_12} / {l21_plus_l23}) * P_1 = {c2_ratio:.4f} * P_1")

    # From (3): lambda_31 * P_3 = lambda_23 * P_2
    # P_3 = (lambda_23 / lambda_31) * P_2
    c3_ratio_p2 = lambda_23 / lambda_31
    # P_3 = c3_ratio_p2 * c2_ratio * P_1
    c3_ratio_p1 = c3_ratio_p2 * c2_ratio
    print(f"From (3): P_3 = ({lambda_23} / {lambda_31}) * P_2 = {c3_ratio_p2:.4f} * P_2 = {c3_ratio_p1:.4f} * P_1")
    print("-" * 30)

    print("Step 3: Substitute into the normalization equation P_0 + P_1 + P_2 + P_3 = 1.")
    # P_0 + P_1 + P_2 + P_3 = (c0_ratio * P_1) + P_1 + (c2_ratio * P_1) + (c3_ratio_p1 * P_1) = 1
    # P_1 * (c0_ratio + 1 + c2_ratio + c3_ratio_p1) = 1
    sum_of_coeffs = c0_ratio + 1 + c2_ratio + c3_ratio_p1
    print(f"( {c0_ratio:.4f} + 1 + {c2_ratio:.4f} + {c3_ratio_p1:.4f} ) * P_1 = 1")
    print(f"{sum_of_coeffs:.4f} * P_1 = 1")
    P1 = 1 / sum_of_coeffs
    print(f"P_1 = {P1:.6f}")
    print("-" * 30)
    
    print("Step 4: Calculate P_0.")
    P0 = c0_ratio * P1
    print(f"P_0 = {c0_ratio:.4f} * P_1 = {c0_ratio:.4f} * {P1:.6f} = {P0:.6f}")
    print("-" * 30)

    print("Step 5: Calculate the final quantity P_0 + P_1.")
    # P_0 + P_1 = c0_ratio * P_1 + P_1 = (c0_ratio + 1) * P_1
    # P_0 + P_1 = (c0_ratio + 1) / (c0_ratio + 1 + c2_ratio + c3_ratio_p1)
    
    # We can write the expression with the numbers
    final_numerator = c0_ratio + 1
    final_denominator = sum_of_coeffs # which is c0_ratio + 3 in this case
    final_result = P0 + P1

    print("The final calculation is P_0 + P_1 = ((c0_ratio) + 1) * P_1")
    print("Which is equal to ((lambda_10/lambda_01) + 1) / ((lambda_10/lambda_01) + 1 + c2_ratio + c3_ratio)")
    print("Let R = lambda_10 / lambda_01:")
    print(f"P_0 + P_1 = (R + 1) / (R + 1 + {c2_ratio} + {c3_ratio_p1}) = (R + 1) / (R + {1 + c2_ratio + c3_ratio_p1})")
    
    print("\n--- Final Equation with Numbers ---")
    print(f"P_0(inf) + P_1(inf) = ( ({lambda_10}/{lambda_01}) + 1 ) / ( ({lambda_10}/{lambda_01}) + 1 + {c2_ratio} + {c3_ratio_p1} )")
    print(f"P_0(inf) + P_1(inf) = ( {c0_ratio} + 1 ) / ( {c0_ratio} + 1 + {c2_ratio} + {c3_ratio_p1} )")
    print(f"P_0(inf) + P_1(inf) = {final_numerator} / {final_denominator}")
    print(f"P_0(inf) + P_1(inf) = {final_result}")


solve_stationary_distribution()