import math
from fractions import Fraction

def solve_betting_problem():
    """
    Calculates and explains the achieved growth rate (W) and the decrease
    in growth rate (Delta W) for a betting scenario with incorrect beliefs.
    """
    # True probabilities
    p = [Fraction(1, 2), Fraction(1, 4), Fraction(1, 8), Fraction(1, 8)]
    # Incorrectly believed probabilities
    q = [Fraction(1, 4), Fraction(1, 2), Fraction(1, 8), Fraction(1, 8)]
    # Bookmaker odds (k-for-1)
    o = [4, 3, 7, 7]

    print("Part 1: Calculate the achieved doubling rate W")
    print("="*50)
    print("You bet according to your incorrect beliefs (q), so your betting fractions are b_i = q_i.")
    print("The achieved growth rate W is calculated using the true probabilities (p):")
    print("W = p1*ln(q1*o1) + p2*ln(q2*o2) + p3*ln(q3*o3) + p4*ln(q4*o4)\n")

    # Build and print the full expression for W
    w_expr_full = []
    for i in range(len(p)):
        w_expr_full.append(f"({p[i]}) * ln(({q[i]}) * {o[i]})")
    print("Substituting the numbers:")
    print(f"W = {' + '.join(w_expr_full)}")

    # Build and print the simplified expression for W
    w_expr_simplified = []
    for i in range(len(p)):
        w_expr_simplified.append(f"({p[i]}) * ln({q[i] * o[i]})")
    print("\nSimplifying the terms inside the logarithm:")
    print(f"W = {' + '.join(w_expr_simplified)}")
    print("W = (1/2) * ln(1) + (1/4) * ln(3/2) + (1/8) * ln(7/8) + (1/8) * ln(7/8)")
    print("\nSince ln(1) = 0, the final expression for W is:")
    print("W = (1/4) * ln(3/2) + (1/4) * ln(7/8)")

    print("\n" + "="*50)
    print("Part 2: Calculate the decrease in doubling rate, Delta W")
    print("="*50)
    print("The decrease in the growth rate is Delta W = W* - W, where W* is the optimal rate.")
    print("This difference is equal to the Kullback-Leibler divergence D_KL(p || q).")
    print("Delta W = p1*ln(p1/q1) + p2*ln(p2/q2) + p3*ln(p3/q3) + p4*ln(p4/q4)\n")

    # Build and print the full expression for Delta W
    dw_expr_full = []
    for i in range(len(p)):
        dw_expr_full.append(f"({p[i]}) * ln(({p[i]})/({q[i]}))")
    print("Substituting the numbers:")
    print(f"Delta W = {' + '.join(dw_expr_full)}")

    # Build and print the simplified expression for Delta W
    dw_expr_simplified = []
    for i in range(len(p)):
        dw_expr_simplified.append(f"({p[i]}) * ln({p[i]/q[i]})")
    print("\nSimplifying the ratios inside the logarithm:")
    print(f"Delta W = {' + '.join(dw_expr_simplified)}")
    print("Delta W = (1/2) * ln(2) + (1/4) * ln(1/2) + (1/8) * ln(1) + (1/8) * ln(1)")
    print("\nUsing ln(1/2) = -ln(2) and ln(1) = 0:")
    print("Delta W = (1/2) * ln(2) - (1/4) * ln(2)")
    print("\nThe final expression for the decrease in doubling rate is:")
    print("Delta W = (1/4) * ln(2)")

solve_betting_problem()