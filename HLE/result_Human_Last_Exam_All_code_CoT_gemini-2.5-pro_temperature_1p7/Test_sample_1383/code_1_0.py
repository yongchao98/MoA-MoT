import math
from fractions import Fraction

def solve_betting_problem():
    """
    Calculates and explains the results of betting with incorrect probabilities.
    """
    # True probabilities
    p = [Fraction(1, 2), Fraction(1, 4), Fraction(1, 8), Fraction(1, 8)]
    # Incorrectly believed probabilities (used for betting fractions b_i)
    q = [Fraction(1, 4), Fraction(1, 2), Fraction(1, 8), Fraction(1, 8)]
    # Bookmaker odds (x-for-1)
    odds_for_1 = [4, 3, 7, 7]
    # Net odds (o_i)
    o = [x - 1 for x in odds_for_1]
    # Payouts (o_i + 1)
    payouts = [x for x in odds_for_1]

    # --- 1. Calculate the achieved growth rate W ---
    print("--- Achieved Growth Rate W (using incorrect probabilities q_i for betting) ---")
    
    w_terms_str = []
    w_simplified_terms_str = []
    
    for i in range(len(p)):
        w_terms_str.append(f"({p[i]}) * log(({q[i]}) * {payouts[i]})")
        # log(1) term will be zero, others are kept
        if q[i] * payouts[i] != 1:
            w_simplified_terms_str.append(f"({p[i]}) * log({q[i] * payouts[i]})")

    w_calc_str = " + ".join(w_terms_str)
    print(f"W = {w_calc_str}")

    w_intermediate_str = " + ".join([f"({p[i]}) * log({q[i] * payouts[i]})" for i in range(len(p))])
    print(f"W = {w_intermediate_str}")
    
    w_final_str = " + ".join(w_simplified_terms_str)
    print(f"Simplified W = {w_final_str}")
    print("\n")

    # --- 2. Calculate the optimal growth rate W* ---
    print("--- Optimal Growth Rate W* (using true probabilities p_i for betting) ---")

    w_star_terms_str = []
    for i in range(len(p)):
        w_star_terms_str.append(f"({p[i]}) * log(({p[i]}) * {payouts[i]})")

    w_star_calc_str = " + ".join(w_star_terms_str)
    print(f"W* = {w_star_calc_str}")
    
    w_star_final_str = " + ".join([f"({p[i]}) * log({p[i] * payouts[i]})" for i in range(len(p))])
    print(f"W* = {w_star_final_str}")
    print("\n")

    # --- 3. Calculate the decrease in growth rate Delta W ---
    print("--- Decrease in Growth Rate Delta W ---")
    print("Delta W = W* - W")
    print(f"Delta W = ( {w_star_final_str} ) - ( {w_final_str} )")

    # The decrease in growth rate is the Kullback-Leibler divergence D_KL(p || q)
    print("\nThis simplifies to the Kullback-Leibler divergence D_KL(p || q) = sum(p_i * log(p_i / q_i))")
    
    kl_terms_str = []
    kl_sum = 0
    for i in range(len(p)):
        if q[i] > 0:
            ratio = p[i] / q[i]
            kl_terms_str.append(f"({p[i]}) * log({ratio})")
            if ratio > 0:
                kl_sum += p[i] * math.log(ratio)

    kl_calc_str = " + ".join(kl_terms_str)
    print(f"Delta W = {kl_calc_str}")

    # Simplify the KL divergence expression by hand for the final result
    # (1/2)*log(2) + (1/4)*log(1/2) = (1/2)*log(2) - (1/4)*log(2) = (1/4)*log(2)
    delta_w_final_str = "(1/4) * log(2)"
    print(f"Simplified Delta W = {delta_w_final_str}")
    print("\n")
    
    # --- Final Answer ---
    # Simplified W = 1/4 * log(3/2) + 1/4 * log(7/8) = 1/4 * log( (3/2)*(7/8) ) = 1/4 * log(21/16)
    final_w_simplified = "(1/4) * log(21/16)"

    # Print final answer in the required format
    print(f"<<<W = {w_final_str}, Delta W = {delta_w_final_str}>>>")


if __name__ == '__main__':
    solve_betting_problem()
