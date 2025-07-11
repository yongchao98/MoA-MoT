import math

def solve_race_growth_rate():
    """
    Calculates and explains the difference between optimal (W*) and actual (W)
    portfolio growth rates based on the Kelly Criterion.
    """
    # Step 1: Define problem parameters
    p_actual = [0.5, 0.25, 0.25]   # Actual probabilities for A, B, C
    p_believed = [0.25, 0.5, 0.25]  # Mistakenly believed probabilities
    payouts = [4, 3, 3]             # Payout ratios
    odds = [p - 1 for p in payouts] # Net odds (b)

    # --- Calculate the optimal strategy and growth rate (W*) ---
    # With actual probabilities, only betting on A has a positive expected value:
    # EV_A = 0.5 * 3 - 0.5 * 1 = 1.0
    # The optimal fraction to bet on A is f*_A = p - (1-p)/b
    f_star_A = p_actual[0] - (1 - p_actual[0]) / odds[0]

    # W* is the expected log wealth growth using the optimal strategy.
    # W* = P_actual(A_wins)*log(1 + f*_A*b_A) + P_actual(A_loses)*log(1 - f*_A)
    term1_W_star = 1 + f_star_A * odds[0]
    term2_W_star = 1 - f_star_A
    W_star = p_actual[0] * math.log(term1_W_star) + (1 - p_actual[0]) * math.log(term2_W_star)

    # --- Calculate the mistaken strategy and the resulting actual growth rate (W) ---
    # With believed probabilities, only betting on B has a positive expected value:
    # EV_B = 0.5 * 2 - 0.5 * 1 = 0.5
    # The mistaken fraction to bet on B is f_B = p - (1-p)/b
    f_mistake_B = p_believed[1] - (1 - p_believed[1]) / odds[1]

    # W is the expected log wealth growth using the mistaken bet, but with actual probabilities.
    # W = P_actual(A)*log(wealth if A) + P_actual(B)*log(wealth if B) + P_actual(C)*log(wealth if C)
    # If A or C wins, the bet on B is lost. If B wins, the bet on B wins.
    term1_W = 1 - f_mistake_B
    term2_W = 1 + f_mistake_B * odds[1]
    W = (p_actual[0] * math.log(term1_W) +
         p_actual[1] * math.log(term2_W) +
         p_actual[2] * math.log(term1_W))

    # --- Calculate the difference and print the result ---
    difference = W_star - W

    print("Step 1: Optimal Growth Rate (W*)")
    print(f"The optimal strategy is to bet a fraction of {f_star_A:.4f} on Competitor A.")
    print(f"The formula for W* is: P(A)*log(1 + f_A*b_A) + (1-P(A))*log(1 - f_A)")
    print(f"W* = {p_actual[0]} * log(1 + {f_star_A:.4f} * {odds[0]}) + {1-p_actual[0]} * log(1 - {f_star_A:.4f})")
    print(f"W* = {p_actual[0]} * log({term1_W_star:.4f}) + {1-p_actual[0]} * log({term2_W_star:.4f})")
    print(f"W* ≈ {W_star:.5f}\n")

    print("Step 2: Actual Growth Rate (W)")
    print(f"The mistaken strategy is to bet a fraction of {f_mistake_B:.4f} on Competitor B.")
    print(f"The formula for W is: P(A)*log(1-f_B) + P(B)*log(1+f_B*b_B) + P(C)*log(1-f_B)")
    print(f"W = {p_actual[0]}*log(1-{f_mistake_B:.4f}) + {p_actual[1]}*log(1+{f_mistake_B:.4f}*{odds[1]}) + {p_actual[2]}*log(1-{f_mistake_B:.4f})")
    print(f"W = {p_actual[0]}*log({term1_W:.4f}) + {p_actual[1]}*log({term2_W:.4f}) + {p_actual[2]}*log({term1_W:.4f})")
    print(f"W ≈ {W:.5f}\n")

    print("Step 3: Final Calculation (W* - W)")
    print(f"W* - W = {W_star:.5f} - ({W:.5f})")
    print(f"W* - W = {difference:.5f}")

    # Return the final numerical answer in the required format
    return f"<<<{difference:.5f}>>>"

# Execute the function and print the final answer
final_answer = solve_race_growth_rate()
print(final_answer)
