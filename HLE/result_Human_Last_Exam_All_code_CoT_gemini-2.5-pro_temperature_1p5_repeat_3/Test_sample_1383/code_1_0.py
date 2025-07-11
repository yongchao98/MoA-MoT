from fractions import Fraction

def solve_growth_rate_problem():
    """
    Solves the bike racing betting problem by calculating growth rates.
    """
    # Step 1: Define true probabilities (p), incorrect probabilities (q), and decimal odds (d)
    p = [Fraction(1, 2), Fraction(1, 4), Fraction(1, 8), Fraction(1, 8)]
    q = [Fraction(1, 4), Fraction(1, 2), Fraction(1, 8), Fraction(1, 8)]
    d = [Fraction(4), Fraction(3), Fraction(7), Fraction(7)]

    # Step 2: Calculate the optimal growth rate W*
    # Find which bets have a positive expectation with true probabilities p_i * d_i > 1
    # p_1 * d_1 = 1/2 * 4 = 2 > 1.  (Bet)
    # p_2 * d_2 = 1/4 * 3 = 3/4 < 1. (Do not bet)
    # p_3 * d_3 = 1/8 * 7 = 7/8 < 1. (Do not bet)
    # p_4 * d_4 = 1/8 * 7 = 7/8 < 1. (Do not bet)
    # The optimal strategy is to only bet on Bike 1.

    p1, d1 = p[0], d[0]
    b_star_1 = (p1 * d1 - 1) / (d1 - 1)  # Kelly fraction: (2-1)/(4-1) = 1/3

    # W* is the expected log growth rate for this strategy.
    # Outcome 1: Bike 1 wins (prob p1 = 1/2).
    # Wealth ratio = (1 - b_star_1) + b_star_1 * d1 = (1 - 1/3) + (1/3)*4 = 2/3 + 4/3 = 2
    # Outcome 2: Bike 1 loses (prob 1-p1 = 1/2).
    # Wealth ratio = 1 - b_star_1 = 1 - 1/3 = 2/3
    # W* = p1 * log(2) + (1 - p1) * log(2/3)
    #    = 1/2 * log(2) + 1/2 * (log(2) - log(3))
    #    = log(2) - 1/2 * log(3)
    w_star_log2_coeff = Fraction(1)
    w_star_log3_coeff = Fraction(-1, 2)

    # Step 3: Calculate the achieved growth rate W with incorrect beliefs
    # User decides bets based on incorrect probabilities q_i * d_i > 1
    # q_1 * d_1 = 1/4 * 4 = 1.      (Do not bet, zero edge)
    # q_2 * d_2 = 1/2 * 3 = 3/2 > 1.  (Bet)
    # q_3 * d_3 = 1/8 * 7 = 7/8 < 1.  (Do not bet)
    # q_4 * d_4 = 1/8 * 7 = 7/8 < 1.  (Do not bet)
    # The user's strategy is to only bet on Bike 2.

    q2, d2 = q[1], d[1]
    b2 = (q2 * d2 - 1) / (d2 - 1) # Kelly fraction: (3/2 - 1)/(3-1) = (1/2)/2 = 1/4

    # W is the expected log growth using user's bet b2, but with true probabilities.
    p2 = p[1] # True probability of Bike 2 winning is 1/4
    # Outcome 1: Bike 2 wins (true prob p2 = 1/4).
    # Wealth ratio = (1 - b2) + b2 * d2 = (1 - 1/4) + (1/4)*3 = 3/4 + 3/4 = 3/2
    # Outcome 2: Bike 2 loses (true prob 1-p2 = 3/4).
    # Wealth ratio = 1 - b2 = 1 - 1/4 = 3/4
    # W = p2 * log(3/2) + (1 - p2) * log(3/4)
    #   = 1/4 * (log(3) - log(2)) + 3/4 * (log(3) - log(4))
    #   = 1/4*log(3) - 1/4*log(2) + 3/4*log(3) - 3/4*log(2^2)
    #   = log(3) - 1/4*log(2) - (3/4)*2*log(2)
    #   = log(3) - 7/4*log(2)
    w_log2_coeff = Fraction(-7, 4)
    w_log3_coeff = Fraction(1)

    # Step 4: Calculate the decrease in growth rate Delta W
    # Delta W = W* - W
    delta_w_log2_coeff = w_star_log2_coeff - w_log2_coeff # 1 - (-7/4) = 11/4
    delta_w_log3_coeff = w_star_log3_coeff - w_log3_coeff # -1/2 - 1 = -3/2
    
    # Step 5: Format and print the results
    print("The doubling rate you will achieve is W:")
    print(f"W = p_2 * log(3/2) + (1 - p_2) * log(3/4)")
    print(f"W = (1/4) * log(3/2) + (3/4) * log(3/4)")
    print(f"W = log({int(w_log3_coeff)}) - ({abs(w_log2_coeff.numerator)}/{w_log2_coeff.denominator}) * log(2)")
    print("\nBy how much has your doubling rate decreased, ΔW = W* - W:")
    print(f"ΔW = W* - W = (log(2) - (1/2)log(3)) - (log(3) - (7/4)log(2))")
    print(f"ΔW = ({delta_w_log2_coeff.numerator}/{delta_w_log2_coeff.denominator}) * log(2) - ({abs(delta_w_log3_coeff.numerator)}/{delta_w_log3_coeff.denominator}) * log(3)")

solve_growth_rate_problem()
print("\n<<<W = log(3) - (7/4)log(2), ΔW = (11/4)log(2) - (3/2)log(3)>>>")