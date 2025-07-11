import math

def calculate_growth_difference():
    """
    Calculates the difference between the optimal (W*) and actual (W) growth rates
    for a three-competitor race given a mistaken belief about probabilities.
    """
    # Step 1: Define probabilities and payout odds
    P_actual = [0.5, 0.25, 0.25]   # Actual probabilities for (A, B, C)
    Q_believed = [0.25, 0.5, 0.25] # Believed probabilities for (A, B, C)
    b_odds = [3.0, 2.0, 2.0]       # Net odds b for (A, B, C)
    o_payouts = [b + 1 for b in b_odds] # Payout ratios o for (A, B, C)

    # Step 2: Determine mistaken bet fractions (f) based on believed probabilities (Q)
    f_mistaken = [0.0, 0.0, 0.0]
    for i in range(3):
        # The Kelly criterion is to bet only if expected payout is greater than 1
        if Q_believed[i] * o_payouts[i] > 1:
            f_mistaken[i] = (Q_believed[i] * o_payouts[i] - 1) / b_odds[i]

    # Step 3: Calculate the actual growth rate (W) using f_mistaken and P_actual
    total_fraction_bet_W = sum(f_mistaken)
    W = 0
    for i in range(3):
        # Wealth factor = 1 - (total bet) + (winnings from this outcome)
        wealth_factor = 1 - total_fraction_bet_W + f_mistaken[i] * b_odds[i]
        W += P_actual[i] * math.log(wealth_factor)

    # Step 4: Determine optimal bet fractions (f_star) based on actual probabilities (P)
    f_star = [0.0, 0.0, 0.0]
    for i in range(3):
        if P_actual[i] * o_payouts[i] > 1:
            f_star[i] = (P_actual[i] * o_payouts[i] - 1) / b_odds[i]

    # Step 5: Calculate the optimal growth rate (W_star) using f_star and P_actual
    total_fraction_bet_W_star = sum(f_star)
    W_star = 0
    for i in range(3):
        wealth_factor = 1 - total_fraction_bet_W_star + f_star[i] * b_odds[i]
        W_star += P_actual[i] * math.log(wealth_factor)
        
    # Step 6: Calculate and print the difference W* - W
    result = W_star - W
    
    print(f"The optimal growth rate (W*) is: {W_star}")
    print(f"The actual growth rate with the mistaken bet (W) is: {W}")
    print("\nThe final equation W* - W is:")
    print(f"{W_star} - ({W}) = {result}")

calculate_growth_difference()