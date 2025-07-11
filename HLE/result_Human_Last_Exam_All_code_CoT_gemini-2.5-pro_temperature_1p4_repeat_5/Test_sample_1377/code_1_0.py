import math

def solve_race_growth_rate():
    """
    Calculates the difference between optimal (W*) and actual (W) 
    growth rates in a betting scenario with mistaken probabilities.
    """
    # --- 1. Define Parameters ---
    # True probabilities for (A, B, C)
    p = {'A': 0.5, 'B': 0.25, 'C': 0.25}
    # Mistaken probabilities for (A, B, C)
    q = {'A': 0.25, 'B': 0.5, 'C': 0.25}
    # Net odds (payout is b+1)
    b = {'A': 3.0, 'B': 2.0, 'C': 2.0}
    competitors = ['A', 'B', 'C']

    # --- 2. Calculate Optimal Strategy and W* ---
    print("--- Calculating Optimal Growth Rate (W*) ---")
    f_star = {'A': 0.0, 'B': 0.0, 'C': 0.0}
    optimal_bet_horse = None
    
    # Check expectation E = p*b - (1-p) for each horse
    for horse in competitors:
        expectation = p[horse] * b[horse] - (1 - p[horse])
        if expectation > 0:
            f_star[horse] = expectation / b[horse]
            optimal_bet_horse = horse

    print(f"Optimal strategy based on true probabilities is to bet on Competitor {optimal_bet_horse}.")
    f_optimal = f_star[optimal_bet_horse]
    print(f"Optimal fraction to bet: {f_optimal:.4f}")
    
    # Calculate W* using the optimal fraction and true probabilities
    p_win_optimal = p[optimal_bet_horse]
    w_star = p_win_optimal * math.log(1 + f_optimal * b[optimal_bet_horse]) + (1 - p_win_optimal) * math.log(1 - f_optimal)
    print(f"W* = p({optimal_bet_horse})*log(1 + f*b) + (1-p({optimal_bet_horse}))*log(1 - f)")
    print(f"W* = {p_win_optimal:.2f}*log(1 + {f_optimal:.4f}*{b[optimal_bet_horse]}) + {1-p_win_optimal:.2f}*log(1 - {f_optimal:.4f})")
    print(f"Optimal growth rate W* = {w_star:.4f}\n")

    # --- 3. Calculate Mistaken Strategy and W ---
    print("--- Calculating Actual Growth Rate (W) ---")
    f_actual = {'A': 0.0, 'B': 0.0, 'C': 0.0}
    actual_bet_horse = None

    # Determine strategy based on the *mistaken* probabilities q
    for horse in competitors:
        expectation = q[horse] * b[horse] - (1 - q[horse])
        if expectation > 0:
            f_actual[horse] = expectation / b[horse]
            actual_bet_horse = horse

    print(f"Mistaken strategy based on perceived probabilities is to bet on Competitor {actual_bet_horse}.")
    f_mistaken = f_actual[actual_bet_horse]
    print(f"Fraction bet based on mistaken beliefs: {f_mistaken:.4f}")

    # Calculate actual growth rate W using the mistaken bet fraction but the *true* probabilities
    p_win_actual = p[actual_bet_horse]
    w_actual = p_win_actual * math.log(1 + f_mistaken * b[actual_bet_horse]) + (1 - p_win_actual) * math.log(1 - f_mistaken)
    print(f"W = p({actual_bet_horse})*log(1 + f*b) + (1-p({actual_bet_horse}))*log(1 - f)")
    print(f"W = {p_win_actual:.2f}*log(1 + {f_mistaken:.4f}*{b[actual_bet_horse]}) + {1-p_win_actual:.2f}*log(1 - {f_mistaken:.4f})")
    print(f"Actual growth rate W = {w_actual:.4f}\n")

    # --- 4. Calculate W* - W ---
    print("--- Calculating the Difference ---")
    difference = w_star - w_actual
    print(f"W* - W = {w_star:.4f} - ({w_actual:.4f})")
    print(f"Final Answer: {difference:.4f}")
    
    # Store the result to be returned in the specified format
    global final_answer_value
    final_answer_value = difference

solve_race_growth_rate()