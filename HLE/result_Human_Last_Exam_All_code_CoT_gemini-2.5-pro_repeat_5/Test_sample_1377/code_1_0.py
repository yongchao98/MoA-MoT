import numpy as np

def calculate_growth_rates():
    """
    Calculates the optimal vs. actual growth rate based on the problem's parameters.
    """
    # Step 1: Define probabilities and payouts
    # True probabilities for A, B, C
    p_true = np.array([1/2, 1/4, 1/4])
    # Believed (mistaken) probabilities for A, B, C
    p_believed = np.array([1/4, 1/2, 1/4])
    # Payout ratios are b:1. A payout of 4:1 means for $1 bet, you get $4 back, which is a $3 profit.
    # So, the net odds 'b' are 3 for A, and 2 for B and C.
    b = np.array([3, 2, 2])

    # Step 2: Calculate betting fractions based on the mistaken belief (f_mistaken)
    f_mistaken = np.zeros(3)
    # A bet is only placed if the expected value is positive. E = p*b - (1-p)
    for i in range(3):
        expected_value = p_believed[i] * b[i] - (1 - p_believed[i])
        if expected_value > 0:
            # Kelly fraction for a single bet: f = p - (1-p)/b
            f_mistaken[i] = p_believed[i] - (1 - p_believed[i]) / b[i]

    # Step 3: Calculate the actual growth rate (W) using f_mistaken and p_true
    wealth_outcomes_W = np.zeros(3)
    total_bet_W = np.sum(f_mistaken)
    for i in range(3):
        # If competitor i wins, wealth is W_0 * (1 - total_bet + f_i * (b_i + 1))
        wealth_outcomes_W[i] = 1 - total_bet_W + f_mistaken[i] * (b[i] + 1)
    
    # The growth rate W is the expected value of the log of wealth, using true probabilities
    W = np.sum(p_true * np.log(wealth_outcomes_W))

    # Step 4: Calculate the optimal betting fractions (f_optimal) based on true probabilities
    f_optimal = np.zeros(3)
    for i in range(3):
        expected_value = p_true[i] * b[i] - (1 - p_true[i])
        if expected_value > 0:
            f_optimal[i] = p_true[i] - (1 - p_true[i]) / b[i]

    # Step 5: Calculate the optimal growth rate (W*)
    wealth_outcomes_W_star = np.zeros(3)
    total_bet_W_star = np.sum(f_optimal)
    for i in range(3):
        wealth_outcomes_W_star[i] = 1 - total_bet_W_star + f_optimal[i] * (b[i] + 1)

    W_star = np.sum(p_true * np.log(wealth_outcomes_W_star))

    # Step 6: Calculate and print the difference W* - W
    difference = W_star - W
    
    # Output the final equation with all the numbers
    print(f"Optimal growth rate W* = {W_star}")
    print(f"Actual growth rate W = {W}")
    print(f"The difference W* - W = {W_star} - ({W}) = {difference}")
    print(f"<<<{difference:.6f}>>>")

calculate_growth_rates()