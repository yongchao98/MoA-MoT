import math

def solve_growth_rate_difference():
    """
    Calculates the difference between the optimal (W*) and actual (W)
    growth rates in a horse race betting scenario.
    """
    # Step 1: Define the parameters of the problem
    p_true = {'A': 1/2, 'B': 1/4, 'C': 1/4}
    p_believed = {'A': 1/4, 'B': 1/2, 'C': 1/4}
    # Net odds (payout is b:1, so a $1 bet yields $b in profit)
    odds = {'A': 4, 'B': 3, 'C': 3}

    # Step 2: Calculate the optimal strategy and growth rate (W*)
    # This strategy is based on the TRUE probabilities.
    # We only bet on outcomes where p * b > 1.
    # Edge A: 0.5 * 4 - 1 = 1.0  (Bet)
    # Edge B: 0.25 * 3 - 1 = -0.25 (Don't bet)
    # Edge C: 0.25 * 3 - 1 = -0.25 (Don't bet)
    # Since only competitor A has a positive edge, we only bet on A.
    
    f_star_A = (p_true['A'] * odds['A'] - 1) / odds['A']

    # Optimal growth rate W* is the expected log return with the optimal bet.
    # Wealth if A wins: 1 + f_star_A * odds['A']
    # Wealth if A loses (B or C wins): 1 - f_star_A
    w_star = (p_true['A'] * math.log(1 + f_star_A * odds['A']) +
              p_true['B'] * math.log(1 - f_star_A) +
              p_true['C'] * math.log(1 - f_star_A))

    # Step 3: Calculate the mistaken strategy and the resulting actual growth rate (W)
    # This strategy is based on the MISTAKEN probabilities.
    # Edge A: 0.25 * 4 - 1 = 0.0   (Don't bet)
    # Edge B: 0.5 * 3 - 1 = 0.5   (Bet)
    # Edge C: 0.25 * 3 - 1 = -0.25 (Don't bet)
    # The mistaken belief leads to betting only on B.
    
    f_prime_B = (p_believed['B'] * odds['B'] - 1) / odds['B']

    # The actual growth rate W results from using the mistaken bet f_prime_B,
    # evaluated with the TRUE probabilities.
    # Wealth if B wins: 1 + f_prime_B * odds['B']
    # Wealth if B loses (A or C wins): 1 - f_prime_B
    w_actual = (p_true['A'] * math.log(1 - f_prime_B) +
                p_true['B'] * math.log(1 + f_prime_B * odds['B']) +
                p_true['C'] * math.log(1 - f_prime_B))

    # Step 4: Calculate the final difference
    difference = w_star - w_actual

    # Print the result in the requested equation format
    print("The final calculation for the difference (W* - W) is:")
    print(f"{w_star:.4f} - ({w_actual:.4f}) = {difference:.4f}")

    return difference

# Execute the function to get the final answer for the '<<<' format.
final_answer = solve_growth_rate_difference()
# The value is returned to be used in the final answer block.
# print(f"\n<<< {final_answer:.4f} >>>") # for my own testing