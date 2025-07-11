import sys
from fractions import Fraction

def calculate_p_limit():
    """
    Calculates the limit for p based on scaling arguments for the integral's L^p norm.
    """

    # Dimensions of the subspaces for linear, quadratic, and cubic coefficients
    d_L = 2  # (a1, a2)
    d_Q = 3  # (a3, a4, a5)
    d_C = 4  # (a6, a7, a8, a9)

    # --- Case 1: Quadratic-dominant region ---
    print("Analyzing the quadratic-dominant region:")
    # Decay exponent of the integral I(a) ~ Lambda^{-k}
    k_Q = 1
    # Scaling exponents for lower-order variables
    # a_L ~ Lambda^alpha_L, a_C ~ Lambda^alpha_C
    alpha_L_Q = Fraction(1, 2)
    alpha_C_Q = Fraction(3, 2)

    # Numerator of the expression for p: sum of weighted dimensions
    # Formula: p <= (alpha_L*d_L + alpha_C*d_C + d_Q) / k
    numerator_Q = alpha_L_Q * d_L + alpha_C_Q * d_C + d_Q
    p_limit_Q = numerator_Q / k_Q

    print(f"Decay exponent k = {k_Q}")
    print(f"Dimensions of subspaces (linear, quadratic, cubic): d_L={d_L}, d_Q={d_Q}, d_C={d_C}")
    print(f"Scaling exponents: alpha_L = {alpha_L_Q}, alpha_C = {alpha_C_Q}")
    print(f"The condition for divergence is p <= (alpha_L*d_L + alpha_C*d_C + d_Q) / k")
    print(f"p <= (({alpha_L_Q})*{d_L} + ({alpha_C_Q})*{d_C} + {d_Q}) / {k_Q}")
    print(f"p <= ({alpha_L_Q * d_L} + {alpha_C_Q * d_C} + {d_Q}) / {k_Q}")
    print(f"p <= {numerator_Q} / {k_Q}")
    print(f"p <= {p_limit_Q}\n")

    # --- Case 2: Cubic-dominant region ---
    print("Analyzing the cubic-dominant region:")
    # Decay exponent
    k_C = Fraction(2, 3)
    # Scaling exponents
    # a_L ~ Lambda^alpha_L, a_Q ~ Lambda^alpha_Q
    alpha_L_C = Fraction(1, 3)
    alpha_Q_C = Fraction(2, 3)

    # Numerator of the expression for p
    # Formula: p <= (alpha_L*d_L + alpha_Q*d_Q + d_C) / k
    numerator_C = alpha_L_C * d_L + alpha_Q_C * d_Q + d_C
    p_limit_C = numerator_C / k_C
    
    print(f"Decay exponent k = {k_C}")
    print(f"Dimensions of subspaces (linear, quadratic, cubic): d_L={d_L}, d_Q={d_Q}, d_C={d_C}")
    print(f"Scaling exponents: alpha_L = {alpha_L_C}, alpha_Q = {alpha_Q_C}")
    print(f"The condition for divergence is p <= (alpha_L*d_L + alpha_Q*d_Q + d_C) / k")
    print(f"p <= (({alpha_L_C})*{d_L} + ({alpha_Q_C})*{d_Q} + {d_C}) / ({k_C})")
    print(f"p <= (({alpha_L_C * d_L}) + ({alpha_Q_C * d_Q}) + {d_C}) / ({k_C})")
    print(f"p <= ({numerator_C}) / ({k_C})")
    print(f"p <= {p_limit_C}\n")

    if p_limit_Q == p_limit_C:
        print(f"Both cases yield the same result: {p_limit_Q}.")
        print("The largest p such that the function I is not in L^p(R^9) is the upper bound of this range.")
        final_p = float(p_limit_Q)
        print(f"Final answer: {final_p}")
        return final_p
    else:
        print("Error: The results from different regions do not match.", file=sys.stderr)
        return None

calculate_p_limit()