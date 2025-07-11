import math

def solve_growth_rate_difference():
    """
    Calculates the difference between the optimal (W*) and actual (W)
    growth rates in a horse race betting scenario.
    """
    # Step 1: Define probabilities and net odds
    # True probabilities
    p = {'A': 1/2, 'B': 1/4, 'C': 1/4}
    # Mistaken probabilities
    q = {'A': 1/4, 'B': 1/2, 'C': 1/4}
    # Net odds (payout-1)
    # A: 4:1 payout -> b=3. B,C: 3:1 payout -> b=2.
    b = {'A': 3, 'B': 2, 'C': 2}

    # --- Step 2: Calculate the optimal growth rate W* ---
    
    # Using true probabilities 'p', check which bets have a positive edge
    # Edge exists if p*(b+1) > 1
    # A: (1/2)*(3+1) = 2.0 > 1 (Positive Edge)
    # B: (1/4)*(2+1) = 0.75 < 1 (Negative Edge)
    # C: (1/4)*(2+1) = 0.75 < 1 (Negative Edge)
    # The optimal strategy is to only bet on A.

    f_star_A = p['A'] - (1 - p['A']) / b['A']
    # If A wins (prob p['A']), wealth is multiplied by (1 + f*b)
    # If A loses (prob 1-p['A']), wealth is multiplied by (1 - f)
    win_factor_star = 1 + f_star_A * b['A']
    loss_factor_star = 1 - f_star_A

    W_star = p['A'] * math.log(win_factor_star) + (1 - p['A']) * math.log(loss_factor_star)
    
    # --- Step 3: Calculate the actual growth rate W ---
    
    # Using mistaken probabilities 'q', check which bets the user would make
    # Edge exists if q*(b+1) > 1
    # A: (1/4)*(3+1) = 1.0 (Zero Edge)
    # B: (1/2)*(2+1) = 1.5 > 1 (Positive Edge)
    # C: (1/4)*(2+1) = 0.75 < 1 (Negative Edge)
    # The mistaken strategy is to only bet on B.
    
    f_B = q['B'] - (1 - q['B']) / b['B']
    
    # Calculate the actual growth 'W' using the mistaken bet f_B but with true probabilities 'p'
    # If A wins (prob p['A']), we lose our bet on B. Wealth factor = 1 - f_B
    # If B wins (prob p['B']), we win our bet on B. Wealth factor = 1 + f_B*b['B']
    # If C wins (prob p['C']), we lose our bet on B. Wealth factor = 1 - f_B

    win_factor_B = 1 + f_B * b['B']
    loss_factor_B = 1 - f_B

    W = p['A'] * math.log(loss_factor_B) + \
        p['B'] * math.log(win_factor_B) + \
        p['C'] * math.log(loss_factor_B)

    # --- Step 4: Calculate the difference ---
    
    difference = W_star - W
    
    print("Optimal betting fraction on A (f*_A): 1/2 - (1/2)/3 = 1/3")
    print("Optimal growth rate (W*): 1/2 * log(1 + 1/3*3) + 1/2 * log(1 - 1/3) = 1/2*log(4/3)")
    print(f"W* = {W_star:.5f}\n")
    
    print("Mistaken betting fraction on B (f_B): 1/2 - (1/2)/2 = 1/4")
    print("Actual growth rate (W): 1/2*log(1-1/4) + 1/4*log(1+1/4*2) + 1/4*log(1-1/4) = 3/4*log(3/4) + 1/4*log(3/2)")
    print(f"W = {W:.5f}\n")
    
    # Print final equation with numbers and the final result
    print("Final Calculation:")
    print(f"W* - W = {W_star:.5f} - ({W:.5f})")
    print(f"W* - W = {difference:.5f}")


solve_growth_rate_difference()
# The exact answer can be represented as 1/4 * log(2048/729)
final_exact_value = 1/4 * math.log(2048/729)
print(f"\nFinal exact value is 1/4 * log(2048/729)")
print(f"<<<{final_exact_value:.5f}>>>")