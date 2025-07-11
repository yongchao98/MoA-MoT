import math

def solve_growth_rate_difference():
    """
    Calculates the difference between the optimal growth rate (W*) and the actual
    growth rate (W) based on mistaken probabilities in a three-competitor race.
    """
    # Step 1: Define parameters
    # True probabilities [A, B, C]
    p_true = [1/2, 1/4, 1/4]
    # Believed probabilities [A, B, C]
    p_believed = [1/4, 1/2, 1/4]
    # Payout odds b:1 [A, B, C]
    b = [4, 3, 3]

    # --- Step 2: Calculate Optimal Strategy and Growth Rate (W*) ---

    # A bet is optimal only if p*(b+1) > 1
    # For true probabilities:
    # A: 0.5 * (4+1) = 2.5 > 1 (Profitable)
    # B: 0.25 * (3+1) = 1.0 (Neutral, no edge)
    # C: 0.25 * (3+1) = 1.0 (Neutral, no edge)
    # The optimal strategy is to only bet on Competitor A.

    # Kelly fraction for a single bet: f = p - q/b = p - (1-p)/b
    f_star_A = p_true[0] - (1 - p_true[0]) / b[0]
    
    # Calculate W* using true probabilities and the optimal fraction for A
    # If A wins, wealth is multiplied by (1 + f*b). If A loses, it's multiplied by (1 - f).
    # The probability of loss is P(B) + P(C) = 1 - P(A).
    ret_A_wins_optimal = 1 + f_star_A * b[0]
    ret_A_loses_optimal = 1 - f_star_A
    
    W_star = p_true[0] * math.log(ret_A_wins_optimal) + (1 - p_true[0]) * math.log(ret_A_loses_optimal)

    # --- Step 3: Determine the Bettor's Mistaken Strategy ---
    
    # Check profitability based on believed probabilities:
    # A: 0.25 * (4+1) = 1.25 > 1 (Believes A is profitable)
    # B: 0.5 * (3+1) = 2.0 > 1 (Believes B is profitable)
    # C: 0.25 * (3+1) = 1.0 (Believes C is neutral)
    # The bettor thinks they should bet on A and B. A simultaneous optimization
    # shows the ideal fractions are invalid (one is negative).
    # Thus, the bettor must choose between betting on A only or B only.

    # Case 1 (Belief): Bet on A only
    f_A_only_believed = p_believed[0] - (1 - p_believed[0]) / b[0]
    growth_A_believed = p_believed[0] * math.log(1 + f_A_only_believed * b[0]) + (1 - p_believed[0]) * math.log(1 - f_A_only_believed)

    # Case 2 (Belief): Bet on B only
    f_B_only_believed = p_believed[1] - (1 - p_believed[1]) / b[1]
    growth_B_believed = p_believed[1] * math.log(1 + f_B_only_believed * b[1]) + (1 - p_believed[1]) * math.log(1 - f_B_only_believed)

    # The bettor chooses the strategy with the higher *believed* growth rate.
    if growth_A_believed > growth_B_believed:
        f_actual_bet = [f_A_only_believed, 0, 0]
    else:
        f_actual_bet = [0, f_B_only_believed, 0]

    # --- Step 4: Calculate the Actual Growth Rate (W) ---

    # The bettor chose to bet f_B = 1/3 on B and nothing on A or C.
    fA, fB, fC = f_actual_bet

    # Calculate actual returns based on this mistaken bet
    ret_if_A_wins = 1 + fA * b[0] - fB - fC  # Lose bet on B
    ret_if_B_wins = 1 - fA + fB * b[1] - fC  # Win bet on B
    ret_if_C_wins = 1 - fA - fB + fC * b[2]  # Lose bet on B

    # Calculate actual growth rate W using TRUE probabilities and the mistaken bet
    W_actual = p_true[0] * math.log(ret_if_A_wins) + \
               p_true[1] * math.log(ret_if_B_wins) + \
               p_true[2] * math.log(ret_if_C_wins)

    # --- Step 5: Compute and Display the Final Result ---
    
    difference = W_star - W_actual

    print("Optimal vs. Actual Growth Rate Calculation\n")
    
    print("Optimal Growth Rate (W*):")
    print(f"The optimal strategy is to bet f_A = {p_true[0]} - (1-{p_true[0]})/{b[0]} = {f_star_A:.4f} on A only.")
    print(f"W* = P(A)*log(1 + f_A*b_A) + (P(B)+P(C))*log(1 - f_A)")
    print(f"W* = {p_true[0]}*log(1 + {f_star_A:.4f}*{b[0]}) + {p_true[1]+p_true[2]}*log(1 - {f_star_A:.4f})")
    print(f"W* = {p_true[0]}*log({ret_A_wins_optimal:.4f}) + {1-p_true[0]}*log({ret_A_loses_optimal:.4f}) = {W_star:.6f}\n")

    print("Actual Growth Rate (W):")
    print(f"The mistaken strategy is to bet f_B = {f_B_only_believed:.4f} on B only.")
    print(f"W = P(A)*log(Return if A wins) + P(B)*log(Return if B wins) + P(C)*log(Return if C wins)")
    print(f"W = {p_true[0]}*log({ret_if_A_wins:.4f}) + {p_true[1]}*log({ret_if_B_wins:.4f}) + {p_true[2]}*log({ret_if_C_wins:.4f}) = {W_actual:.6f}\n")

    print("Final Result (W* - W):")
    print(f"The difference in growth rates is:")
    print(f"W* - W = {W_star:.6f} - ({W_actual:.6f}) = {difference:.6f}")
    
    return difference

# Run the calculation and print the final answer in the required format
final_answer = solve_growth_rate_difference()
print(f"\n<<<{final_answer:.6f}>>>")
