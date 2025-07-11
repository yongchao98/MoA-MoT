import math

def calculate_growth_rate_difference():
    """
    Calculates the difference between the optimal (W*) and actual (W) growth rates
    for a three-competitor race with mistaken probabilities.
    """
    # Step 1: Define problem parameters
    # True probabilities
    p = {'A': 1/2, 'B': 1/4, 'C': 1/4}
    # Believed (mistaken) probabilities
    q = {'A': 1/4, 'B': 1/2, 'C': 1/4}
    # Net odds (payout is b:1)
    b = {'A': 3, 'B': 2, 'C': 2}

    # Step 2: Calculate the optimal strategy and W* based on true probabilities (p)
    # Check for positive edge: p*(b+1) - 1 > 0
    edge_A_opt = p['A'] * (b['A'] + 1) - 1
    
    # The optimal strategy is to only bet on A, as it's the only one with a positive edge.
    f_star_A = edge_A_opt / b['A']
    
    # W* = E[ln(Wealth)] = p_A*ln(1 + f*_A*b_A) + (p_B+p_C)*ln(1 - f*_A)
    w_star = p['A'] * math.log(1 + f_star_A * b['A']) + (p['B'] + p['C']) * math.log(1 - f_star_A)

    # Step 3: Calculate the actual strategy and W based on mistaken probabilities (q)
    # Check for positive edge: q*(b+1) - 1 > 0
    edge_B_actual = q['B'] * (b['B'] + 1) - 1
    
    # The mistaken strategy is to only bet on B.
    f_actual_B = edge_B_actual / b['B']
    
    # The actual growth W is calculated using the mistaken bet (f_actual_B) but with true probabilities (p).
    # W = p_B*ln(1 + f_B*b_B) + (p_A+p_C)*ln(1-f_B)
    w_actual = p['B'] * math.log(1 + f_actual_B * b['B']) + (p['A'] + p['C']) * math.log(1 - f_actual_B)
    
    # Step 4: Compute the difference W* - W
    difference = w_star - w_actual

    # --- Output the results ---
    print("The difference W* - W is calculated as follows:\n")

    print("1. Optimal Growth Rate (W*):")
    print("   - Based on true probabilities (p_A=1/2, p_B=1/4, p_C=1/4).")
    print(f"   - Optimal fraction to bet on A: f*_A = ({p['A']} * ({b['A']}+1) - 1) / {b['A']} = {f_star_A:.4f}")
    print("   - Equation: W* = p_A*ln(1 + f*_A*b_A) + (p_B+p_C)*ln(1 - f*_A)")
    print(f"   - W* = ({p['A']})*ln(1 + {f_star_A:.4f}*{b['A']}) + ({p['B']}+{p['C']})*ln(1 - {f_star_A:.4f})")
    print(f"   - W* = {w_star:.5f}\n")

    print("2. Actual Growth Rate (W):")
    print("   - Based on mistaken probabilities (q_A=1/4, q_B=1/2, q_C=1/4).")
    print(f"   - Actual fraction bet on B: f_B = ({q['B']} * ({b['B']}+1) - 1) / {b['B']} = {f_actual_B:.4f}")
    print("   - Equation: W = p_B*ln(1 + f_B*b_B) + (p_A+p_C)*ln(1 - f_B)")
    print(f"   - W = ({p['B']})*ln(1 + {f_actual_B:.4f}*{b['B']}) + ({p['A']}+{p['C']})*ln(1 - {f_actual_B:.4f})")
    print(f"   - W = {w_actual:.5f}\n")

    print("3. Final Result (W* - W):")
    print(f"   W* - W = {w_star:.5f} - ({w_actual:.5f})")
    print(f"   W* - W = {difference:.5f}")

calculate_growth_rate_difference()