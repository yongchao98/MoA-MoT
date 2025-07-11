import math

def solve_race_growth_rate():
    """
    Calculates the difference between optimal and actual growth rates in a betting scenario.
    """

    # Step 1: Define Probabilities and Payouts
    p_true = [0.5, 0.25, 0.25]  # True probabilities for (A, B, C)
    p_believed = [0.25, 0.5, 0.25] # Mistakenly believed probabilities
    # Payout ratios N:1 mean net winnings are N for a 1 unit bet.
    b = [4, 3, 3]  # Net payout odds for (A, B, C)

    print("--- Problem Setup ---")
    print(f"True probabilities P(A,B,C): ({p_true[0]}, {p_true[1]}, {p_true[2]})")
    print(f"Believed probabilities P'(A,B,C): ({p_believed[0]}, {p_believed[1]}, {p_believed[2]})")
    print(f"Payout odds b(A,B,C): ({b[0]}:1, {b[1]}:1, {b[2]}:1)\n")

    # Step 2: Calculate the Optimal Strategy (f*) and Growth Rate (W*)
    # We use the TRUE probabilities. A bet is optimal if p*(b+1) > 1.
    # For A: 0.5 * (4 + 1) = 2.5 > 1 (Bet)
    # For B: 0.25 * (3 + 1) = 1.0 (Do not bet)
    # For C: 0.25 * (3 + 1) = 1.0 (Do not bet)
    # The optimal strategy is to only bet on A.

    f_star_A = p_true[0] - (1 - p_true[0]) / b[0]
    
    print("--- Calculating Optimal Growth Rate W* ---")
    print("Based on true probabilities, only betting on A is optimal.")
    print(f"Optimal fraction to bet on A, f_A* = P(A) - (1 - P(A)) / b_A")
    print(f"f_A* = {p_true[0]} - (1 - {p_true[0]}) / {b[0]} = {f_star_A:.4f}")

    # W* is the growth rate using f_star_A and true probabilities.
    # W* = P(A)*log(1 + f_A* * b_A) + (1-P(A))*log(1 - f_A*)
    term1_star = 1 + f_star_A * b[0]
    term2_star = 1 - f_star_A
    W_star = p_true[0] * math.log(term1_star) + (1 - p_true[0]) * math.log(term2_star)

    print("\nThe optimal growth rate W* is calculated as:")
    print("W* = P(A) * log(1 + f_A* * b_A) + (1-P(A)) * log(1 - f_A*)")
    print(f"W* = {p_true[0]} * log(1 + {f_star_A:.4f} * {b[0]}) + {1-p_true[0]} * log(1 - {f_star_A:.4f})")
    print(f"W* = {p_true[0]} * log({term1_star:.4f}) + {1-p_true[0]} * log({term2_star:.4f})")
    print(f"W* = {W_star:.4f}\n")

    # Step 3: Calculate the Mistaken Strategy (f) and Actual Growth Rate (W)
    # We use the BELIEVED probabilities. A bet is placed if p'*(b+1) > 1.
    # For A: 0.25 * (4 + 1) = 1.25 > 1 (Bet)
    # For B: 0.5 * (3 + 1) = 2.0 > 1 (Bet)
    # For C: 0.25 * (3 + 1) = 1.0 (Do not bet)
    # The mistaken strategy is to bet on A and B. The fractions were solved analytically.
    
    f_A = 7/44
    f_B = 17/44

    print("--- Calculating Actual Growth Rate W ---")
    print("Based on mistaken probabilities, betting on A and B seems optimal.")
    print(f"The calculated fractions are f_A = 7/44 = {f_A:.4f} and f_B = 17/44 = {f_B:.4f}")
    
    # W is the growth rate using mistaken fractions f_A, f_B and TRUE probabilities.
    # W = P(A)*log(1-f_B+f_A*b_A) + P(B)*log(1-f_A+f_B*b_B) + P(C)*log(1-f_A-f_B)
    term_A_wins = 1 - f_B + f_A * b[0]
    term_B_wins = 1 - f_A + f_B * b[1]
    term_C_wins = 1 - f_A - f_B
    W = p_true[0] * math.log(term_A_wins) + p_true[1] * math.log(term_B_wins) + p_true[2] * math.log(term_C_wins)

    print("\nThe actual growth rate W is calculated using these fractions but with true probabilities:")
    print("W = P(A)*log(1-f_B+f_A*b_A) + P(B)*log(1-f_A+f_B*b_B) + P(C)*log(1-f_A-f_B)")
    print(f"W = {p_true[0]}*log(1-{f_B:.4f}+{f_A:.4f}*{b[0]}) + {p_true[1]}*log(1-{f_A:.4f}+{f_B:.4f}*{b[1]}) + {p_true[2]}*log(1-{f_A:.4f}-{f_B:.4f})")
    print(f"W = {p_true[0]}*log({term_A_wins:.4f}) + {p_true[1]}*log({term_B_wins:.4f}) + {p_true[2]}*log({term_C_wins:.4f})")
    print(f"W = {W:.4f}\n")

    # Step 4: Compute the Difference
    result = W_star - W
    
    print("--- Final Calculation ---")
    print("The difference is W* - W")
    print(f"W* - W = {W_star:.4f} - ({W:.4f})")
    print(f"W* - W = {result:.4f}")
    
    # Final answer in the required format
    print(f"\n<<<{result:.4f}>>>")

solve_race_growth_rate()