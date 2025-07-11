import math
from fractions import Fraction

def solve_betting_problem():
    """
    Calculates and prints the doubling rate W and its decrease ΔW
    based on the given probabilities and odds.
    """
    # True probabilities
    p = [Fraction(1, 2), Fraction(1, 4), Fraction(1, 8), Fraction(1, 8)]
    
    # Incorrectly believed probabilities (used as betting fractions b)
    q = [Fraction(1, 4), Fraction(1, 2), Fraction(1, 8), Fraction(1, 8)]
    
    # Bookmaker odds (k-for-1, so payoff is k)
    o = [4, 3, 7, 7]
    
    # --- Calculate Achieved Doubling Rate W ---
    # W = Σ p_i * log(q_i * o_i)
    
    print("The achieved doubling rate W is calculated as:")
    w_expr = "W = "
    w_terms_str = []
    for i in range(len(p)):
        term_val = q[i] * o[i]
        w_terms_str.append(f"({p[i]}) * ln({q[i]} * {o[i]})")
    w_expr += " + ".join(w_terms_str)
    print(w_expr)

    w_simplified_expr = "W = "
    w_simplified_terms_str = []
    for i in range(len(p)):
        term_val = q[i] * o[i]
        if term_val != 1: # log(1) is 0, so we can omit it
             w_simplified_terms_str.append(f"({p[i]}) * ln({term_val})")
    w_simplified_expr += " + ".join(w_simplified_terms_str)
    print("\nSimplifying the expression for W:")
    print(w_simplified_expr)
    print("-" * 30)

    # --- Calculate the Decrease in Doubling Rate ΔW ---
    # ΔW = W* - W = D_KL(p || q) = Σ p_i * log(p_i / q_i)
    
    print("The decrease in doubling rate ΔW is calculated as:")
    delta_w_expr = "ΔW = "
    delta_w_terms_str = []
    for i in range(len(p)):
        ratio = p[i] / q[i]
        delta_w_terms_str.append(f"({p[i]}) * ln(({p[i]})/({q[i]}))")
    delta_w_expr += " + ".join(delta_w_terms_str)
    print(delta_w_expr)
    
    delta_w_simplified_expr = "ΔW = "
    delta_w_value = 0
    
    # Calculate the sum for the final compact expression
    p1, p2, _, _ = p
    q1, q2, _, _ = q
    # ΔW = p1*log(p1/q1) + p2*log(p2/q2) = (1/2)log(2) + (1/4)log(1/2) = (1/4)log(2)
    final_coeff = p[0] * math.log(p[0]/q[0]) + p[1] * math.log(p[1]/q[1])
    # Compare to (1/4)log(2)
    final_coeff_fraction = Fraction(final_coeff / math.log(2)).limit_denominator()

    print("\nSimplifying the expression for ΔW:")
    print(f"ΔW = ({p[0]}) * ln({p[0]/q[0]}) + ({p[1]}) * ln({p[1]/q[1]})")
    print(f"ΔW = ({final_coeff_fraction}) * ln(2)")

# Run the calculation and print the results
solve_betting_problem()

# Final Answer
# W = (1/4) * ln(3/2) + (1/4) * ln(7/8)
# ΔW = (1/4) * ln(2)
final_answer_str = "W = (1/4) * ln(3/2) + (1/4) * ln(7/8), ΔW = (1/4) * ln(2)"
print(f"\n<<<{final_answer_str}>>>")