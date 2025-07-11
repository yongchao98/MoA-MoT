import math

def calculate_growth_rate_difference():
    """
    Calculates the difference between the optimal Kelly growth rate (W*) and
    the actual growth rate (W) from betting based on mistaken probabilities.
    """
    # Step 0: Define the problem parameters
    # Competitors: A, B, C
    # p_true: True probabilities of winning
    # p_believed: Mistakenly believed probabilities of winning
    # b_odds: Net odds for a winning bet (payout - stake) / stake
    p_true = [0.5, 0.25, 0.25]
    p_believed = [0.25, 0.5, 0.25]
    # Payout ratio X:1 means net odds are X-1. So 4:1 -> b=3, 3:1 -> b=2
    b_odds = [3, 2, 2]

    # --- Step 1: Calculate the optimal strategy and growth rate (W*) ---
    # The optimal strategy is based on the TRUE probabilities.
    # The Kelly fraction f = p - q/b, where q = 1-p. Bet only if f > 0.
    f_star = [0, 0, 0]
    # Check Competitor A
    if p_true[0] * (b_odds[0] + 1) > 1:
        f_star[0] = p_true[0] - (1 - p_true[0]) / b_odds[0]
    # In this case, only A has a positive edge: 0.5 - (0.5)/3 = 1/3

    f_star_total = sum(f_star)
    
    # Calculate returns for each outcome under optimal betting
    R_star_A = 1 - f_star_total + f_star[0] * (b_odds[0] + 1)
    R_star_B = 1 - f_star_total + f_star[1] * (b_odds[1] + 1)
    R_star_C = 1 - f_star_total + f_star[2] * (b_odds[2] + 1)

    # W* is the expected log-return, using true probabilities
    W_star = p_true[0] * math.log(R_star_A) + \
             p_true[1] * math.log(R_star_B) + \
             p_true[2] * math.log(R_star_C)

    # --- Step 2: Calculate the mistaken strategy and resulting growth rate (W) ---
    # The mistaken strategy is based on the BELIEVED probabilities.
    f_actual = [0, 0, 0]
    # Check Competitor B using believed probability
    if p_believed[1] * (b_odds[1] + 1) > 1:
        f_actual[1] = p_believed[1] - (1 - p_believed[1]) / b_odds[1]
    # In this case, only B has a positive edge from the mistaken view: 0.5 - 0.5/2 = 1/4

    f_actual_total = sum(f_actual)

    # Calculate actual returns for each outcome under mistaken betting
    R_actual_A = 1 - f_actual_total + f_actual[0] * (b_odds[0] + 1)
    R_actual_B = 1 - f_actual_total + f_actual[1] * (b_odds[1] + 1)
    R_actual_C = 1 - f_actual_total + f_actual[2] * (b_odds[2] + 1)

    # W is the expected log-return of the mistaken strategy, evaluated with TRUE probabilities
    W_actual = p_true[0] * math.log(R_actual_A) + \
               p_true[1] * math.log(R_actual_B) + \
               p_true[2] * math.log(R_actual_C)

    # --- Step 3: Calculate the difference W* - W and print results ---
    result = W_star - W_actual

    print(f"Optimal Growth Rate (W*):")
    print(f"W* = P_true(A)*ln(R*_A) + P_true(B)*ln(R*_B) + P_true(C)*ln(R*_C)")
    print(f"W* = {p_true[0]} * ln({R_star_A:.4f}) + {p_true[1]} * ln({R_star_B:.4f}) + {p_true[2]} * ln({R_star_C:.4f}) = {W_star:.4f}")
    
    print("\nActual Growth Rate (W):")
    print(f"W = P_true(A)*ln(R_A) + P_true(B)*ln(R_B) + P_true(C)*ln(R_C)")
    print(f"W = {p_true[0]} * ln({R_actual_A:.4f}) + {p_true[1]} * ln({R_actual_B:.4f}) + {p_true[2]} * ln({R_actual_C:.4f}) = {W_actual:.4f}")

    # The difference W* - W can be expressed in simplified analytical form
    # W* = ln(2) - 0.5*ln(3)
    # W = ln(3) - 1.75*ln(2)
    # W* - W = 2.75*ln(2) - 1.5*ln(3)
    term1_val = 2.75 * math.log(2)
    term2_val = 1.5 * math.log(3)

    print("\nThe difference W* - W can be simplified to the equation:")
    print("W* - W = 2.75 * ln(2) - 1.5 * ln(3)")
    print(f"W* - W = {term1_val:.4f} - {term2_val:.4f}")
    print(f"The final difference is: {result:.4f}")
    
# Execute the calculation
if __name__ == "__main__":
    calculate_growth_rate_difference()