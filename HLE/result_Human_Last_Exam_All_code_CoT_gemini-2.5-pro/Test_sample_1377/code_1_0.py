import math

def solve_kelly_problem():
    """
    Calculates the difference between the optimal (W*) and actual (W) 
    Kelly growth rates based on a mistaken belief in probabilities.
    """
    # --- Step 1: Define Probabilities and Payouts ---
    # True probabilities
    p_true = {'A': 1/2, 'B': 1/4, 'C': 1/4}
    # Mistaken probabilities
    p_mistaken = {'A': 1/4, 'B': 1/2, 'C': 1/4}
    # Net odds (b). Payout "x:1" means net odds b = x-1.
    # A is 4:1 -> b_A = 3
    # B and C are 3:1 -> b_B = b_C = 2
    b = {'A': 3, 'B': 2, 'C': 2}

    print("This problem calculates the difference between the optimal growth rate (W*) and the actual growth rate (W) from a mistaken betting strategy.")
    print("-" * 80)

    # --- Step 2 & 3: Calculate Optimal Strategy and Growth Rate (W*) ---
    print("\n1. Calculating the Optimal Strategy and Growth Rate (W*)...")
    print(f"   True Probabilities P = (A: {p_true['A']:.2f}, B: {p_true['B']:.2f}, C: {p_true['C']:.2f})")
    print(f"   Net Odds b = (A: {b['A']}, B: {b['B']}, C: {b['C']})")

    # The optimal strategy is to bet only on A, as it's the only favorable bet
    # p_A(b_A+1) = 0.5 * (3+1) = 2.0 > 1
    # p_B(b_B+1) = 0.25 * (2+1) = 0.75 < 1
    f_A_star = (p_true['A'] * (b['A'] + 1) - 1) / b['A']
    
    print(f"\n   The only favorable bet is on A. Optimal fraction f_A* = {f_A_star:.4f}")

    # Calculate W*
    W_star = (p_true['A'] * math.log(1 + b['A'] * f_A_star) +
              p_true['B'] * math.log(1 - f_A_star) +
              p_true['C'] * math.log(1 - f_A_star))

    print("   The optimal growth rate W* is calculated as:")
    print(f"   W* = P(A)*log(1 + b_A*f_A*) + P(B)*log(1 - f_A*) + P(C)*log(1 - f_A*)")
    print(f"   W* = {p_true['A']:.2f}*log(1 + {b['A']}*{f_A_star:.4f}) + {p_true['B']:.2f}*log(1 - {f_A_star:.4f}) + {p_true['C']:.2f}*log(1 - {f_A_star:.4f})")
    print(f"   W* = {W_star:.4f}")
    print("-" * 80)
    
    # --- Step 4 & 5: Calculate Mistaken Strategy and Actual Growth Rate (W) ---
    print("\n2. Calculating the Mistaken Strategy and Actual Growth Rate (W)...")
    print(f"   Mistaken Probabilities P' = (A: {p_mistaken['A']:.2f}, B: {p_mistaken['B']:.2f}, C: {p_mistaken['C']:.2f})")
    
    # The mistaken strategy is to bet only on B
    # p'_A(b_A+1) = 0.25 * (3+1) = 1.0 (Not strictly positive edge)
    # p'_B(b_B+1) = 0.5 * (2+1) = 1.5 > 1
    f_B_actual = (p_mistaken['B'] * (b['B'] + 1) - 1) / b['B']

    print(f"\n   Based on P', the only favorable bet is on B. Fraction f_B = {f_B_actual:.4f}")
    
    # Calculate W (actual) using true probabilities and mistaken fractions
    W_actual = (p_true['A'] * math.log(1 - f_B_actual) +
                p_true['B'] * math.log(1 + b['B'] * f_B_actual) +
                p_true['C'] * math.log(1 - f_B_actual))
    
    print("   The actual growth rate W is calculated using P (true probs) and f_B (mistaken bet):")
    print(f"   W = P(A)*log(1 - f_B) + P(B)*log(1 + b_B*f_B) + P(C)*log(1 - f_B)")
    print(f"   W = {p_true['A']:.2f}*log(1 - {f_B_actual:.4f}) + {p_true['B']:.2f}*log(1 + {b['B']}*{f_B_actual:.4f}) + {p_true['C']:.2f}*log(1 - {f_B_actual:.4f})")
    print(f"   W = {W_actual:.4f}")
    print("-" * 80)
    
    # --- Step 6: Compute and Print the Difference ---
    difference = W_star - W_actual
    
    print("\n3. Final Calculation: W* - W")
    print("   The difference in growth rates is:")
    print(f"   W* - W = {W_star:.4f} - ({W_actual:.4f})")
    print(f"   W* - W = {difference:.4f}")

if __name__ == '__main__':
    solve_kelly_problem()
<<<0.2582>>>