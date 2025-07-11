import math

def solve_growth_rate_difference():
    """
    Calculates the difference between the optimal growth rate (W*) and the actual
    growth rate (W) from a mistaken betting strategy in a three-competitor race.
    """

    # --- Step 1: Define Probabilities and Payouts ---
    # True probabilities
    p_true = {'A': 1/2, 'B': 1/4, 'C': 1/4}
    # Net payout ratios (b:1)
    b = {'A': 4, 'B': 3, 'C': 3}

    # --- Step 2: Calculate the Optimal Growth Rate (W*) ---
    # Based on maximizing E[log(W)] with true probabilities, the optimal
    # strategy is to bet f_A = 3/8 on A, and f_B = f_C = 0.
    f_star_A = 3/8
    
    # Wealth factor if A wins
    w_star_A = 1 + f_star_A * b['A']
    # Wealth factor if A loses (B or C wins)
    w_star_loss = 1 - f_star_A
    
    # W* = P(A)log(W_A) + P(B)log(W_loss) + P(C)log(W_loss)
    # W* = P(A)log(W_A) + (P(B)+P(C))log(W_loss)
    w_star = p_true['A'] * math.log(w_star_A) + (p_true['B'] + p_true['C']) * math.log(w_star_loss)
    
    # An equivalent simpler form is W* = 0.5 * log(25/16) = log(5/4)

    # --- Step 3: Calculate the Actual Growth Rate (W) ---
    # Fractions derived from maximizing E[log(W)] with mistaken probabilities
    # P' = (1/4, 1/2, 1/4). This yields f'_A = 7/44, f'_B = 17/44.
    f_mistaken_A = 7/44
    f_mistaken_B = 17/44

    # Calculate the wealth factors for each outcome using the mistaken fractions
    # Wealth if A wins
    w_actual_A = 1 - f_mistaken_B + f_mistaken_A * b['A']
    # Wealth if B wins
    w_actual_B = 1 - f_mistaken_A + f_mistaken_B * b['B']
    # Wealth if C wins (no bet placed, so just lose the stakes on A and B)
    w_actual_C = 1 - f_mistaken_A - f_mistaken_B
    
    # Calculate the actual growth rate W using TRUE probabilities and mistaken fractions
    w_actual = (p_true['A'] * math.log(w_actual_A) +
                p_true['B'] * math.log(w_actual_B) +
                p_true['C'] * math.log(w_actual_C))

    # --- Step 4: Compute and Print the Difference ---
    difference = w_star - w_actual

    # As requested, print the final equation with numerical values
    # W* - W = log(5/4) - [1/2 * log(W_A_actual) + 1/4 * log(W_B_actual) + 1/4 * log(W_C_actual)]
    # W* - W = log(1.25) - (0.5 * log(1.25) + 0.25 * log(2.0) + 0.25 * log(0.4545...))
    # W* - W = 0.5 * log(1.25) - 0.25 * log(2.0) - 0.25 * log(0.4545...)
    
    term1_val = 0.5 * math.log(w_star_A / w_star_loss) # This is 0.5*log(1.25)
    term2_val = p_true['B'] * math.log(w_actual_B)
    term3_val = p_true['C'] * math.log(w_actual_C)

    # Note: w_star - w_actual can be simplified to:
    # 0.5*log(1.25) - 0.25*log(2) - 0.25*log(5/11)
    
    print("Equation for the difference:")
    print(f"W* - W = 0.5 * log({w_star_A / w_star_loss}) - 0.25 * log({w_actual_B}) - 0.25 * log({w_actual_C})")
    print(f"W* - W = {0.5 * math.log(1.25):.4f} - {0.25 * math.log(2.0):.4f} - {0.25 * math.log(w_actual_C):.4f}")
    print(f"W* - W = {difference:.4f}")

solve_growth_rate_difference()
print(f"\n<<<{(math.log(5/4) - (0.5 * math.log(5/4) + 0.25 * math.log(2) + 0.25 * math.log(5/11))):.4f}>>>")
