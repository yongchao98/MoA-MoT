import math
from fractions import Fraction

def solve_betting_problem():
    """
    Solves the Kelly criterion betting problem step-by-step, calculating W*, W, and ΔW.
    The results are printed as equations with fractions and natural logs.
    """
    # --- Problem Setup ---
    # True probabilities
    p = {
        1: Fraction(1, 2),
        2: Fraction(1, 4),
        3: Fraction(1, 8),
        4: Fraction(1, 8)
    }
    # Incorrect (believed) probabilities
    q = {
        1: Fraction(1, 4),
        2: Fraction(1, 2),
        3: Fraction(1, 8),
        4: Fraction(1, 8)
    }
    # Bookmaker odds (d-for-1 payoff)
    d = {
        1: 4,
        2: 3,
        3: 7,
        4: 7
    }

    print("This script calculates the doubling rates for a bike race betting scenario.")
    print("-" * 50)

    # --- Part 1: Optimal Doubling Rate (W*) ---
    print("Part 1: Calculating the Optimal Doubling Rate (W*)\n")
    print("We use the true probabilities (p_i) to find the optimal betting fractions (b_i*).")
    print("A bet is placed on bike 'i' only if the expected value is positive (p_i * d_i > 1).")

    b_star = {i: Fraction(0) for i in range(1, 5)}
    for i in sorted(p.keys()):
        val = p[i] * d[i]
        print(f"For Bike {i}: p_{i} * d_{i} = ({p[i]}) * {d[i]} = {val}")
        if val > 1:
            print(f"  -> Favorable bet. Calculating Kelly fraction b_{i}*.")
            b_star[i] = (p[i] * d[i] - 1) / (d[i] - 1)
            print(f"     b_{i}* = (p_{i}*d_{i} - 1) / (d_{i} - 1) = ({val} - 1) / ({d[i]} - 1) = {b_star[i]}")
        else:
            print(f"  -> Not a favorable bet. b_{i}* = 0.")

    total_bet_star = sum(b_star.values())
    print(f"\nThe optimal strategy is to bet a fraction of {b_star[1]} on Bike 1 and 0 on the others.")
    
    print("\nThe optimal doubling rate W* is W* = SUM[ p_i * log(1 - SUM(b_j*) + b_i* * d_i) ]")
    
    log_arg_win = 1 - total_bet_star + b_star[1] * d[1]
    log_arg_lose = 1 - total_bet_star

    print(f"\nBreaking it down by outcome:")
    print(f"If Bike 1 wins (prob {p[1]}): log argument is 1 - {total_bet_star} + {b_star[1]}*{d[1]} = {log_arg_win}")
    print(f"If any other bike wins (prob {1-p[1]}): log argument is 1 - {total_bet_star} + 0 = {log_arg_lose}")
    
    print(f"\nSo, the equation for W* is:")
    print(f"W* = ({p[1]})*log({log_arg_win}) + ({p[2]})*log({log_arg_lose}) + ({p[3]})*log({log_arg_lose}) + ({p[4]})*log({log_arg_lose})")
    print(f"W* = (1/2)*log(2) + (1/2)*log(2/3)")
    print(f"   (which simplifies to log(2) - (1/2)*log(3))")

    print("\n" + "-" * 50)

    # --- Part 2: Achieved Doubling Rate (W) ---
    print("Part 2: Calculating the Achieved Doubling Rate (W) with Incorrect Beliefs\n")
    print("The bettor uses incorrect probabilities (q_i) to determine their betting fractions (b_i).")
    print("A bet is placed on bike 'i' only if the bettor believes q_i * d_i > 1.")

    b_incorrect = {i: Fraction(0) for i in range(1, 5)}
    for i in sorted(q.keys()):
        val = q[i] * d[i]
        print(f"For Bike {i}: q_{i} * d_{i} = ({q[i]}) * {d[i]} = {val}")
        if val > 1:
            print(f"  -> Bettor believes this is a favorable bet. Calculating b_{i}.")
            b_incorrect[i] = (q[i] * d[i] - 1) / (d[i] - 1)
            print(f"     b_{i} = (q_{i}*d_{i} - 1) / (d_{i} - 1) = ({val} - 1) / ({d[i]} - 1) = {b_incorrect[i]}")
        else:
            print(f"  -> Bettor believes this is not a favorable bet. b_{i} = 0.")
    
    total_bet_incorrect = sum(b_incorrect.values())
    print(f"\nThe bettor's strategy is to bet a fraction of {b_incorrect[2]} on Bike 2 and 0 on the others.")
    
    print("\nThe achieved doubling rate W uses true probabilities (p_i) but the bettor's fractions (b_i):")
    print("W = SUM[ p_i * log(1 - SUM(b_j) + b_i * d_i) ]")

    log_arg_win_b2 = 1 - total_bet_incorrect + b_incorrect[2] * d[2]
    log_arg_lose_b2 = 1 - total_bet_incorrect

    print(f"\nBreaking it down by outcome (using true probabilities):")
    print(f"If Bike 2 wins (prob {p[2]}): log argument is 1 - {total_bet_incorrect} + {b_incorrect[2]}*{d[2]} = {log_arg_win_b2}")
    print(f"If any other bike wins: log argument is 1 - {total_bet_incorrect} + 0 = {log_arg_lose_b2}")

    print(f"\nSo, the equation for W is:")
    print(f"W = ({p[1]})*log({log_arg_lose_b2}) + ({p[2]})*log({log_arg_win_b2}) + ({p[3]})*log({log_arg_lose_b2}) + ({p[4]})*log({log_arg_lose_b2})")
    print(f"W = (3/4)*log(3/4) + (1/4)*log(3/2)")
    print(f"   (which simplifies to log(3) - (7/4)*log(2))")

    print("\n" + "-" * 50)

    # --- Part 3: Decrease in Doubling Rate (ΔW) ---
    print("Part 3: Calculating the Decrease in Doubling Rate (ΔW)\n")
    print("ΔW = W* - W")
    print("W* = log(2) - (1/2)*log(3)")
    print("W  = log(3) - (7/4)*log(2)")
    print("\nΔW = [log(2) - (1/2)*log(3)] - [log(3) - (7/4)*log(2)]")
    print("ΔW = (1 + 7/4)*log(2) + (-1/2 - 1)*log(3)")
    print("ΔW = (11/4)*log(2) - (3/2)*log(3)")

solve_betting_problem()