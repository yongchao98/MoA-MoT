import math

def calculate_growth_difference():
    """
    Calculates the difference between the optimal (W*) and actual (W) wealth growth rates
    based on a mistaken belief in race probabilities.
    """
    
    # --- Problem Parameters ---
    # True probabilities of winning
    p_true = {'A': 1/2, 'B': 1/4, 'C': 1/4}
    # Net payout odds (e.g., 4:1 means odds are 4)
    odds = {'A': 4, 'B': 3, 'C': 3}

    # --- Step 1: Calculate the Actual Growth Rate (W) from the Mistaken Bet ---

    # Based on a mistaken belief of probabilities (1/4, 1/2, 1/4), the calculated
    # optimal fractions to bet would be f_A = 7/44 and f_B = 17/44.
    f_mistaken = {'A': 7/44, 'B': 17/44, 'C': 0}

    # Calculate the wealth multipliers for each possible outcome using these mistaken betting fractions.
    # If A wins: your wealth is multiplied by (1 + 4*f_A - f_B - f_C)
    # If B wins: your wealth is multiplied by (1 + 3*f_B - f_A - f_C)
    # If C wins: your wealth is multiplied by (1 - f_A - f_B) since f_C=0
    w_outcome_A_mistaken = 1 + odds['A'] * f_mistaken['A'] - f_mistaken['B'] - f_mistaken['C']
    w_outcome_B_mistaken = 1 + odds['B'] * f_mistaken['B'] - f_mistaken['A'] - f_mistaken['C']
    w_outcome_C_mistaken = 1 + odds['C'] * f_mistaken['C'] - f_mistaken['A'] - f_mistaken['B']
    
    # Calculate the actual growth rate 'W' using the TRUE probabilities.
    W = (p_true['A'] * math.log(w_outcome_A_mistaken) +
         p_true['B'] * math.log(w_outcome_B_mistaken) +
         p_true['C'] * math.log(w_outcome_C_mistaken))


    # --- Step 2: Calculate the Optimal Growth Rate (W*) ---

    # The truly optimal strategy, based on the true probabilities, is to bet only on
    # Competitor A. The optimal fraction is f_A = p_true_A - (1-p_true_A)/odds_A.
    f_optimal_A = p_true['A'] - (1 - p_true['A']) / odds['A'] # 0.5 - (0.5 / 4) = 3/8

    # Calculate wealth multipliers for the truly optimal bet.
    # If A wins: your wealth is multiplied by (1 + 4*f_A)
    # If B or C wins: your wealth is multiplied by (1 - f_A)
    w_outcome_A_optimal = 1 + odds['A'] * f_optimal_A
    w_outcome_not_A_optimal = 1 - f_optimal_A

    # Calculate the optimal growth rate 'W*' using the TRUE probabilities.
    # Note: P(Not A) = P(B) + P(C)
    W_star = (p_true['A'] * math.log(w_outcome_A_optimal) +
              (p_true['B'] + p_true['C']) * math.log(w_outcome_not_A_optimal))


    # --- Step 3: Display the final result ---
    difference = W_star - W

    print("This script calculates the difference between the optimal and actual wealth growth rates.\n")
    print(f"The calculated optimal growth rate (W*) is: {W_star:.4f}")
    print(f"The actual growth rate (W) with the mistaken bet is: {W:.4f}\n")
    print("The difference (W* - W) is therefore:")
    print(f"{W_star:.4f} - {W:.4f} = {difference:.4f}")

calculate_growth_difference()