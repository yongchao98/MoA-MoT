import math
from fractions import Fraction
from collections import defaultdict

def solve_betting_problem():
    """
    Calculates and prints the results for the betting problem.
    """
    # Given probabilities and odds
    p_true = [Fraction(1, 2), Fraction(1, 4), Fraction(1, 8), Fraction(1, 8)]
    p_believed = [Fraction(1, 4), Fraction(1, 2), Fraction(1, 8), Fraction(1, 8)]
    odds = [4, 3, 7, 7]
    num_outcomes = len(p_true)

    # --- 1. Calculate the optimal growth rate W* ---
    # Optimal bets b_star are the true probabilities p_true
    bets_optimal = p_true
    
    w_star_terms = defaultdict(Fraction)
    for i in range(num_outcomes):
        payoff = bets_optimal[i] * odds[i]
        w_star_terms[payoff] += p_true[i]

    w_star_parts = []
    for payoff, coeff in sorted(w_star_terms.items()):
        if coeff > 0:
            w_star_parts.append(f"({coeff}) * ln({payoff})")
    w_star_expression = " + ".join(w_star_parts)

    print("The optimal expected logarithmic growth rate W* is given by:")
    print(f"W* = Sum(p_i * ln(p_i * o_i))")
    print(f"W* = {w_star_expression}\n")

    # --- 2. Calculate the achieved growth rate W ---
    # Bets b are the incorrectly believed probabilities p_believed
    bets_incorrect = p_believed

    w_terms = defaultdict(Fraction)
    for i in range(num_outcomes):
        payoff = bets_incorrect[i] * odds[i]
        # We use the true probabilities p_true for the expectation
        w_terms[payoff] += p_true[i]
        
    w_parts = []
    for payoff, coeff in sorted(w_terms.items()):
        # ln(1) = 0, so we can omit this term for clarity
        if coeff > 0 and payoff != 1:
            w_parts.append(f"({coeff}) * ln({payoff})")
    w_expression = " + ".join(w_parts)

    print("The growth rate W achieved with the incorrect probabilities is:")
    print(f"W = Sum(p_i * ln(q_i * o_i))")
    print(f"W = {w_expression}\n")

    # --- 3. Calculate the decrease in growth rate Delta W ---
    # Delta_W = W* - W = D_KL(p_true || p_believed)
    
    delta_w_terms = defaultdict(Fraction)
    total_delta_w_numeric = 0
    
    for i in range(num_outcomes):
        ratio = p_true[i] / p_believed[i]
        delta_w_terms[ratio] += p_true[i]
        if ratio > 0:
            total_delta_w_numeric += p_true[i] * math.log(ratio)

    delta_w_parts = []
    for ratio, coeff in sorted(delta_w_terms.items()):
        # ln(1) = 0, so omit it
        if coeff > 0 and ratio != 1:
             delta_w_parts.append(f"({coeff}) * ln({ratio})")
    
    # The sum simplifies nicely
    final_delta_w_coeff = Fraction(0)
    for ratio, coeff in delta_w_terms.items():
        if ratio == 2:
            final_delta_w_coeff += coeff
        elif ratio == Fraction(1,2):
            final_delta_w_coeff -= coeff
            
    delta_w_expression_sum = " + ".join(delta_w_parts)
    delta_w_expression_simplified = f"({final_delta_w_coeff}) * ln(2)"

    print("The decrease in the growth rate, ΔW = W* - W, is:")
    print(f"ΔW = Sum(p_i * ln(p_i / q_i))")
    print(f"   = {delta_w_expression_sum}")
    print(f"   = {delta_w_expression_simplified}")

    # Prepare final answer for grading format
    final_w_star = w_star_expression
    final_w = w_expression
    final_delta_w = delta_w_expression_simplified
    
    # Printing the final combined answer as requested.
    # print(f"\n<<<W* = {final_w_star}, W = {final_w}, ΔW = {final_delta_w}>>>")

solve_betting_problem()