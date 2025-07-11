import numpy as np

def solve_growth_rates():
    """
    Calculates and compares optimal vs. actual growth rates based on the Kelly criterion.
    """
    # --- 1. Define problem parameters ---
    # True probabilities for (A, B, C)
    p_true = np.array([1/2, 1/4, 1/4])
    # Mistakenly believed probabilities
    p_believed = np.array([1/4, 1/2, 1/4])
    # Payout odds (b-to-1, so b is net winnings for a 1 unit bet)
    b = np.array([4, 3, 3])

    # --- 2. Calculate betting fractions based on MISTAKEN beliefs (f_mistaken) ---
    # Based on mistaken beliefs, the edge for C is zero, so we only bet on A and B.
    # The maximization of the Kelly formula G = sum(p_i * log(wealth_factor_i))
    # leads to a system of linear equations for the betting fractions f_A and f_B:
    # 37*f_A - 23*f_B = -3
    # 17*f_A + 37*f_B = 17
    A_matrix = np.array([[37, -23], [17, 37]])
    B_vector = np.array([-3, 17])
    f_A_mistaken, f_B_mistaken = np.linalg.solve(A_matrix, B_vector)
    f_mistaken = np.array([f_A_mistaken, f_B_mistaken, 0])

    # --- 3. Calculate the ACTUAL growth rate (W) ---
    # This uses the mistaken fractions but the TRUE probabilities of outcomes.
    # The wealth factor if outcome `i` wins is: 1 + b_i*f_i - sum(f_j for j!=i)
    growth_factors_mistaken = np.array([
        1 + b[0] * f_mistaken[0] - f_mistaken[1] - f_mistaken[2],  # A wins
        1 - f_mistaken[0] + b[1] * f_mistaken[1] - f_mistaken[2],  # B wins
        1 - f_mistaken[0] - f_mistaken[1] + b[2] * f_mistaken[2]   # C wins
    ])
    # The actual expected log growth rate (W) is the sum of p_true * log(growth_factor)
    W = np.sum(p_true * np.log(growth_factors_mistaken))

    # --- 4. Calculate the OPTIMAL betting fractions (f_optimal) based on TRUE probabilities ---
    # Check the edge for each bet with true probabilities: p*b - (1-p)
    # Only competitor A has a positive edge (1/2 * 4 - 1/2 = 1.5 > 0).
    # The others have zero edge (1/4 * 3 - 3/4 = 0). So we only bet on A.
    f_A_optimal = p_true[0] - (1 - p_true[0]) / b[0]
    f_optimal = np.array([f_A_optimal, 0, 0])

    # --- 5. Calculate the OPTIMAL growth rate (W*) ---
    # This uses the optimal fractions and true probabilities.
    growth_factors_optimal = np.array([
        1 + b[0] * f_optimal[0] - f_optimal[1] - f_optimal[2],  # A wins
        1 - f_optimal[0] + b[1] * f_optimal[1] - f_optimal[2],  # B wins
        1 - f_optimal[0] - f_optimal[1] + b[2] * f_optimal[2]   # C wins
    ])
    # The optimal expected log growth rate (W*)
    W_star = np.sum(p_true * np.log(growth_factors_optimal))

    # --- 6. Calculate and print the final result ---
    difference = W_star - W
    
    # Per the instructions, outputting each number in the final equation.
    print(f"{W_star:.6f} - {W:.6f} = {difference:.6f}")

solve_growth_rates()