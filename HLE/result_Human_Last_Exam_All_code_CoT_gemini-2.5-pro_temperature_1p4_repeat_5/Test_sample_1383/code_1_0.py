import math
from fractions import Fraction

def solve_betting_problem():
    """
    Calculates and explains the doubling rates for a bike race betting scenario.
    """
    # Step 1: Define the given probabilities and odds.
    # p are the true probabilities, q are the incorrectly believed probabilities.
    p = [Fraction(1, 2), Fraction(1, 4), Fraction(1, 8), Fraction(1, 8)]
    q = [Fraction(1, 4), Fraction(1, 2), Fraction(1, 8), Fraction(1, 8)]
    # o are the 'for-1' odds, which correspond to decimal odds.
    o = [4, 3, 7, 7]

    print("This script calculates the achieved doubling rate (W) and its decrease from the optimal rate (ΔW).")
    print("All calculations are in terms of natural logs (log) and fractions.\n")

    # Step 2: Calculate the achieved doubling rate W.
    # The user bets based on incorrect probabilities q, so the betting fractions b_i = q_i.
    # The actual expected growth W is calculated using the true probabilities p.
    # Formula: W = Σ p_i * log(q_i * o_i)
    print("--- 1. Calculating the Achieved Doubling Rate W ---")
    print("You bet according to your incorrect beliefs (q), so your betting fractions are:")
    print(f"b = [{', '.join(map(str, q))}]")
    print("The achieved rate W is the expectation over the true probabilities (p):")
    
    w_eq_parts = []
    w_simp_parts = []
    for i in range(len(p)):
        w_eq_parts.append(f"({p[i]}) * log(({q[i]}) * {o[i]})")
        term_val = q[i] * o[i]
        w_simp_parts.append(f"({p[i]}) * log({term_val})")
        
    print(f"W = {' + '.join(w_eq_parts)}")
    print(f"W = {' + '.join(w_simp_parts)}")
    print("Simplifying the expression for W:")
    # W = (1/2)*log(1) + (1/4)*log(3/2) + (1/8)*log(7/8) + (1/8)*log(7/8)
    # Since log(1) = 0 and we can combine the last two terms:
    print("W = (1/4) * log(3/2) + (1/4) * log(7/8)\n")

    # Step 3: Calculate the optimal doubling rate W*.
    # This is achieved by betting with the true probabilities p.
    # Formula: W* = Σ p_i * log(p_i * o_i)
    print("--- 2. Calculating the Optimal Doubling Rate W* ---")
    print("The optimal strategy is to bet according to the true probabilities (p).")
    
    w_star_simp_parts = []
    for i in range(len(p)):
        term_val = p[i] * o[i]
        w_star_simp_parts.append(f"({p[i]}) * log({term_val})")
        
    print(f"W* = {' + '.join(w_star_simp_parts)}")
    print("Simplifying the expression for W*:")
    # W* = (1/2)*log(2) + (1/4)*log(3/4) + (1/8)*log(7/8) + (1/8)*log(7/8)
    print("W* = (1/2) * log(2) + (1/4) * log(3/4) + (1/4) * log(7/8)\n")
    
    # Step 4: Calculate the decrease in the doubling rate, ΔW.
    # ΔW = W* - W
    print("--- 3. Calculating the Decrease in Doubling Rate ΔW ---")
    print("ΔW is the difference W* - W. This is also the KL-divergence D_KL(p || q).")
    print("ΔW = [ (1/2)log(2) + (1/4)log(3/4) + (1/4)log(7/8) ] - [ (1/4)log(3/2) + (1/4)log(7/8) ]")
    print("ΔW = (1/2)log(2) + (1/4)log(3/4) - (1/4)log(3/2)")
    print("Using log properties, this simplifies to the KL-divergence, Σ p_i * log(p_i / q_i):")

    delta_w_parts = []
    for i in range(len(p)):
        ratio = p[i] / q[i]
        delta_w_parts.append(f"({p[i]}) * log({ratio})")
    
    print(f"ΔW = {' + '.join(delta_w_parts)}")
    # ΔW = (1/2)*log(2) + (1/4)*log(1/2) + (1/8)*log(1) + (1/8)*log(1)
    # ΔW = (1/2)*log(2) - (1/4)*log(2) + 0 + 0
    print("ΔW = (1/2) * log(2) - (1/4) * log(2)")
    print("ΔW = (1/4) * log(2)")

if __name__ == '__main__':
    solve_betting_problem()