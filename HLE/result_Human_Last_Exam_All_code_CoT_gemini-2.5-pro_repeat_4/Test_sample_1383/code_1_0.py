from fractions import Fraction

def solve_betting_problem():
    """
    Calculates and explains the optimal and achieved logarithmic growth rates for a betting scenario.
    """
    # p: true probabilities
    p = [Fraction(1, 2), Fraction(1, 4), Fraction(1, 8), Fraction(1, 8)]
    # q: incorrect (believed) probabilities
    q = [Fraction(1, 4), Fraction(1, 2), Fraction(1, 8), Fraction(1, 8)]
    # d: decimal odds (payoff multipliers, e.g., 4-for-1 means d=4)
    d = [Fraction(4), Fraction(3), Fraction(7), Fraction(7)]

    print("--- 1. Optimal Doubling Rate (W*) ---")
    print("This is the rate achieved by betting according to the true probabilities p.")
    print("Formula: W* = Σ p_i * ln(p_i * d_i)")

    # Build and print the equation for W* with all numbers
    w_star_terms = []
    w_star_values = []
    for i in range(len(p)):
        w_star_terms.append(f"({p[i]}) * ln(({p[i]}) * {d[i]})")
        term_val = p[i] * d[i]
        w_star_values.append(f"({p[i]}) * ln({term_val})")

    print(f"\nW* = {' + '.join(w_star_terms)}")
    print(f"W* = {' + '.join(w_star_values)}")
    # Simplified expression for W*
    w_star_final_str = "(1/2) * ln(2) + (1/4) * ln(3/4) + (1/4) * ln(7/8)"
    print(f"Simplified W* = {w_star_final_str}")

    print("\n" + "="*50 + "\n")

    print("--- 2. Achieved Doubling Rate (W) ---")
    print("This is the rate achieved by betting according to the incorrect probabilities q.")
    print("Formula: W = Σ p_i * ln(q_i * d_i)")

    # Build and print the equation for W with all numbers
    w_terms = []
    w_values = []
    for i in range(len(p)):
        w_terms.append(f"({p[i]}) * ln(({q[i]}) * {d[i]})")
        term_val = q[i] * d[i]
        w_values.append(f"({p[i]}) * ln({term_val})")

    print(f"\nW = {' + '.join(w_terms)}")
    print(f"W = {' + '.join(w_values)}")
    
    # Simplified expression for W
    w_final_str = "(1/4) * ln(3/2) + (1/4) * ln(7/8)"
    print(f"Simplified W = {w_final_str}")


    print("\n" + "="*50 + "\n")

    print("--- 3. Decrease in Doubling Rate (ΔW) ---")
    print("This is the difference W* - W, also known as the KL Divergence D_KL(p || q).")
    print("Formula: ΔW = Σ p_i * ln(p_i / q_i)")

    # Build and print the equation for ΔW with all numbers
    delta_w_terms = []
    delta_w_values = []
    for i in range(len(p)):
        ratio = p[i] / q[i]
        delta_w_terms.append(f"({p[i]}) * ln(({p[i]})/({q[i]}))")
        delta_w_values.append(f"({p[i]}) * ln({ratio})")

    print(f"\nΔW = {' + '.join(delta_w_terms)}")
    print(f"ΔW = {' + '.join(delta_w_values)}")
    
    # Further simplification of ΔW
    print("ΔW = (1/2)*ln(2) + (1/4)*ln(1/2) + (1/8)*ln(1) + (1/8)*ln(1)")
    print("ΔW = (1/2)*ln(2) - (1/4)*ln(2)")
    delta_w_final_str = "(1/4) * ln(2)"
    print(f"Final ΔW = {delta_w_final_str}")

    # Final answer in the required format
    final_answer = f"W = {w_final_str}, ΔW = {delta_w_final_str}"
    print("\n" + "="*50 + "\n")
    print(f"<<<{final_answer}>>>")


solve_betting_problem()