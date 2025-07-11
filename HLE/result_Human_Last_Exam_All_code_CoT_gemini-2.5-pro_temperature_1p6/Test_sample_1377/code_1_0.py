import math

def solve_growth_rate_difference():
    """
    Calculates and prints the difference between optimal (W*) and actual (W)
    wealth growth rates based on the Kelly Criterion.
    """
    # Step 1: Define problem parameters
    p_true = {'A': 1/2, 'B': 1/4, 'C': 1/4}
    p_believed = {'A': 1/4, 'B': 1/2, 'C': 1/4}
    
    # Payout ratios are 4:1, 3:1, 3:1.
    # Net odds b = payout - 1
    odds = {'A': 3, 'B': 2, 'C': 2}

    # Step 2: Calculate the optimal strategy and growth rate (W*)
    # The Kelly criterion suggests betting only if p * (b+1) > 1.
    # For the true probabilities:
    # A: (1/2) * (3 + 1) = 2 > 1. Bet.
    # B: (1/4) * (2 + 1) = 0.75 < 1. Do not bet.
    # C: (1/4) * (2 + 1) = 0.75 < 1. Do not bet.
    # So, the optimal strategy is to bet only on A.

    p_A = p_true['A']
    b_A = odds['A']
    # Optimal fraction for a single bet: f = p - (1-p)/b
    f_A_star = p_A - (1 - p_A) / b_A

    # Calculate W* using this optimal fraction under true probabilities
    # W* = p_A * log(1 + b_A*f) + (1-p_A) * log(1 - f)
    # The (1-p_A) term accounts for B or C winning.
    w_star = p_A * math.log(1 + b_A * f_A_star) + (1 - p_A) * math.log(1 - f_A_star)

    # Step 3: Calculate the actual strategy and growth rate (W)
    # The mistaken strategy is based on the believed probabilities:
    # A: (1/4) * (3 + 1) = 1. Do not bet (zero edge).
    # B: (1/2) * (2 + 1) = 1.5 > 1. Bet.
    # C: (1/4) * (2 + 1) = 0.75 < 1. Do not bet.
    # So, the mistaken strategy is to bet only on B.

    p_B_believed = p_believed['B']
    b_B = odds['B']
    # Fraction bet on B is calculated with believed probability
    f_B_actual = p_B_believed - (1 - p_B_believed) / b_B

    # Calculate the actual growth rate W using the mistaken bet f_B_actual
    # but under the true outcome probabilities.
    # If A wins (prob p_true['A']), wealth is (1 - f_B_actual).
    # If B wins (prob p_true['B']), wealth is (1 + b_B * f_B_actual).
    # If C wins (prob p_true['C']), wealth is (1 - f_B_actual).
    
    w_actual = (
        p_true['A'] * math.log(1 - f_B_actual) +
        p_true['B'] * math.log(1 + b_B * f_B_actual) +
        p_true['C'] * math.log(1 - f_B_actual)
    )

    # Step 4: Compute the difference and print the results
    
    print("This script calculates the difference between the optimal and actual wealth growth rates.\n")
    
    print("--- Optimal Growth Rate (W*) ---")
    print("The optimal strategy is to bet only on competitor A.")
    print(f"Optimal fraction to bet on A (f_A*): {p_A:.2f} - (1 - {p_A:.2f}) / {b_A} = {f_A_star:.4f}")
    # W* = (1/2)*ln(2) + (1/2)*ln(2/3) = ln(2) - (1/2)*ln(3)
    print(f"W* = p_A*ln(1 + b_A*f_A*) + (1-p_A)*ln(1-f_A*)")
    print(f"W* = {p_A:.2f}*ln(1 + {b_A}*{f_A_star:.4f}) + {1-p_A:.2f}*ln(1 - {f_A_star:.4f})")
    print(f"W* = (1/2)*ln(2) - (1/2)*ln(3) = {w_star:.4f}\n")
    
    print("--- Actual Growth Rate (W) ---")
    print("The mistaken strategy is to bet only on competitor B.")
    print(f"Mistaken fraction bet on B (f_B): {p_B_believed:.2f} - (1 - {p_B_believed:.2f}) / {b_B} = {f_B_actual:.4f}")
    print("W = p_A*ln(1-f_B) + p_B*ln(1+b_B*f_B) + p_C*ln(1-f_B)")
    print(f"W = {p_true['A']:.2f}*ln(1-{f_B_actual:.4f}) + {p_true['B']:.2f}*ln(1+{b_B}*{f_B_actual:.4f}) + {p_true['C']:.2f}*ln(1-{f_B_actual:.4f})")
    # W = (3/4)*ln(3/4) + (1/4)*ln(3/2) = ln(3) - (7/4)*ln(2)
    print(f"W = ln(3) - (7/4)*ln(2) = {w_actual:.4f}\n")
    
    print("--- Difference (W* - W) ---")
    difference = w_star - w_actual
    # Final symbolic expression: W* - W = (ln(2) - (1/2)ln(3)) - (ln(3) - (7/4)ln(2)) = (11/4)ln(2) - (3/2)ln(3)
    print("W* - W = (ln(2) - (1/2)ln(3)) - (ln(3) - (7/4)ln(2))")
    print("W* - W = (11/4)ln(2) - (3/2)ln(3)")
    print(f"W* - W = {w_star:.4f} - ({w_actual:.4f}) = {difference:.4f}")
    
    return difference

# Execute the function and capture the final answer for the user.
final_answer = solve_growth_rate_difference()
# The 'answer' format is specified in the prompt for internal processing.
# The user sees the detailed printout from the function above.