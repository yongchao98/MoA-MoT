from fractions import Fraction

def solve_growth_rate_problem():
    """
    Calculates and explains the achieved and optimal growth rates, and their difference.
    """
    # Step 1: Define the given probabilities and odds.
    p_true = [Fraction(1, 2), Fraction(1, 4), Fraction(1, 8), Fraction(1, 8)]
    q_believed = [Fraction(1, 4), Fraction(1, 2), Fraction(1, 8), Fraction(1, 8)]
    # "4-for-1" odds means you get 4 back for a 1 bet, so decimal odds are 4, 3, 7, 7.
    d_odds = [4, 3, 7, 7]

    print("--- Calculating the Achieved Doubling Rate (W) ---")
    # The bettor uses their incorrect beliefs (q) as their betting fractions (b).
    b_mistaken = q_believed
    # The actual growth rate W is calculated using true probabilities (p) and the mistaken bets (b).
    # W = p1*ln(b1*d1) + p2*ln(b2*d2) + p3*ln(b3*d3) + p4*ln(b4*d4)
    w_eq_parts = []
    for i in range(len(p_true)):
        w_eq_parts.append(f"({p_true[i]}) * ln({b_mistaken[i]} * {d_odds[i]})")
    
    print("The formula for W is:")
    print(f"W = {' + '.join(w_eq_parts)}")

    # Calculate the values inside the logarithms
    w_calc_parts = []
    for i in range(len(p_true)):
        w_calc_parts.append(f"({p_true[i]}) * ln({b_mistaken[i] * d_odds[i]})")
    
    print("\nPlugging in the numbers:")
    print(f"W = {' + '.join(w_calc_parts)}")
    
    print("\nSince ln(1) = 0, the first term disappears. The final expression for W is:")
    final_W = "(1/4) * ln(3/2) + (1/4) * ln(7/8)"
    print(f"W = {final_W}\n")

    print("--- Calculating the Optimal Doubling Rate (W*) ---")
    # The optimal betting fractions (b*) are the true probabilities (p).
    b_optimal = p_true
    # W* = p1*ln(p1*d1) + p2*ln(p2*d2) + p3*ln(p3*d3) + p4*ln(p4*d4)
    w_star_eq_parts = []
    for i in range(len(p_true)):
        w_star_eq_parts.append(f"({p_true[i]}) * ln({b_optimal[i]} * {d_odds[i]})")
    
    print("The formula for the optimal rate W* is:")
    print(f"W* = {' + '.join(w_star_eq_parts)}")

    # Calculate the values inside the logarithms
    w_star_calc_parts = []
    for i in range(len(p_true)):
        w_star_calc_parts.append(f"({p_true[i]}) * ln({b_optimal[i] * d_odds[i]})")
    
    print("\nPlugging in the numbers:")
    print(f"W* = {' + '.join(w_star_calc_parts)}")

    print("\nThe final expression for W* is:")
    final_W_star = "(1/2) * ln(2) + (1/4) * ln(3/4) + (1/4) * ln(7/8)"
    print(f"W* = {final_W_star}\n")

    print("--- Calculating the Decrease in Doubling Rate (ΔW) ---")
    print("The decrease is the difference ΔW = W* - W:")
    print(f"ΔW = [ {final_W_star} ] - [ {final_W} ]")
    print("Simplifying the expression by canceling terms:")
    print("ΔW = (1/2) * ln(2) + (1/4) * ln(3/4) - (1/4) * ln(3/2)")
    print("Using log properties ln(a/b) = ln(a) - ln(b) and ln(4) = 2*ln(2):")
    print("ΔW = (1/2)*ln(2) + (1/4)*(ln(3) - ln(4)) - (1/4)*(ln(3) - ln(2))")
    print("ΔW = (1/2)*ln(2) + (1/4)ln(3) - (1/4)*(2*ln(2)) - (1/4)ln(3) + (1/4)ln(2)")
    print("ΔW = (1/2)*ln(2) - (1/2)*ln(2) + (1/4)*ln(2)")
    print("The final result for the decrease is:")
    final_delta_W = "(1/4) * ln(2)"
    print(f"ΔW = {final_delta_W}")

    # Final answer in the required format
    print(f"\n<<<W = {final_W}, ΔW = {final_delta_W}>>>")

solve_growth_rate_problem()