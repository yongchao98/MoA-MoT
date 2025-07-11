import math

def solve_kelly_problem():
    """
    Calculates the difference between optimal and actual wealth growth rates
    based on a misunderstanding of race probabilities.
    """
    # Step 1: Define Probabilities and Payouts
    # True probabilities
    p_A, p_B, p_C = 1/2, 1/4, 1/4
    # Mistaken probabilities
    q_A, q_B, q_C = 1/4, 1/2, 1/4
    # Payouts (decimal odds o) and net odds (b)
    o_A, o_B, o_C = 4, 3, 3
    b_A, b_B, b_C = o_A - 1, o_B - 1, o_C - 1

    print("--- Calculating Optimal Growth Rate (W*) ---")
    # Step 2: Determine optimal strategy using true probabilities (p)
    # A bet is only optimal if p * o > 1
    # For A: p_A * o_A = 0.5 * 4 = 2.0 > 1 (Bet)
    # For B: p_B * o_B = 0.25 * 3 = 0.75 < 1 (Don't bet)
    # For C: p_C * o_C = 0.25 * 3 = 0.75 < 1 (Don't bet)
    # The optimal strategy is to only bet on Competitor A.

    # Calculate the optimal Kelly fraction for A
    f_A_star = p_A - (1 - p_A) / b_A
    print(f"Optimal strategy: Bet only on A.")
    print(f"Optimal fraction f_A* = P(A) - (1 - P(A)) / b_A = {p_A} - (1 - {p_A}) / {b_A} = {f_A_star:.4f}")

    # Calculate the optimal growth rate W*
    # W* = p_A * ln(1 + f_A* * b_A) + (p_B + p_C) * ln(1 - f_A*)
    w_star_val1 = math.log(1 + f_A_star * b_A)
    w_star_val2 = math.log(1 - f_A_star)
    W_star = p_A * w_star_val1 + (p_B + p_C) * w_star_val2
    print("\nOptimal growth rate equation (W*):")
    print(f"W* = P(A)*ln(1 + f_A* * b_A) + (P(B) + P(C))*ln(1 - f_A*)")
    print(f"W* = {p_A}*ln(1 + {f_A_star:.4f} * {b_A}) + ({p_B} + {p_C})*ln(1 - {f_A_star:.4f})")
    print(f"W* = {p_A}*ln({1 + f_A_star * b_A:.4f}) + {p_B + p_C}*ln({1 - f_A_star:.4f}) = {W_star:.4f}")
    
    print("\n" + "="*50 + "\n")

    print("--- Calculating Actual Growth Rate (W) ---")
    # Step 3: Determine mistaken strategy using mistaken probabilities (q)
    # A bet is made if q * o > 1
    # For A: q_A * o_A = 0.25 * 4 = 1.0 (Don't bet)
    # For B: q_B * o_B = 0.5 * 3 = 1.5 > 1 (Bet)
    # For C: q_C * o_C = 0.25 * 3 = 0.75 < 1 (Don't bet)
    # The mistaken strategy is to only bet on Competitor B.

    # Calculate the mistaken Kelly fraction for B
    f_B_prime = q_B - (1 - q_B) / b_B
    print(f"Mistaken strategy: Bet only on B.")
    print(f"Mistaken fraction f_B' = Q(B) - (1 - Q(B)) / b_B = {q_B} - (1 - {q_B}) / {b_B} = {f_B_prime:.4f}")

    # Calculate the actual growth rate W using true probabilities (p) and mistaken fraction (f_B')
    # W = p_B * ln(1 + f_B' * b_B) + (p_A + p_C) * ln(1 - f_B')
    w_val1 = math.log(1 + f_B_prime * b_B)
    w_val2 = math.log(1 - f_B_prime)
    W = p_B * w_val1 + (p_A + p_C) * w_val2
    print("\nActual growth rate equation (W):")
    print(f"W = P(B)*ln(1 + f_B' * b_B) + (P(A) + P(C))*ln(1 - f_B')")
    print(f"W = {p_B}*ln(1 + {f_B_prime:.4f} * {b_B}) + ({p_A} + {p_C})*ln(1 - {f_B_prime:.4f})")
    print(f"W = {p_B}*ln({1 + f_B_prime * b_B:.4f}) + {p_A + p_C}*ln({1 - f_B_prime:.4f}) = {W:.4f}")

    print("\n" + "="*50 + "\n")

    # Step 4: Calculate the difference W* - W
    difference = W_star - W
    print("--- Final Result ---")
    print(f"Difference = W* - W")
    print(f"Difference = {W_star:.4f} - ({W:.4f})")
    print(f"Difference = {difference:.4f}")

    return difference

if __name__ == '__main__':
    result = solve_kelly_problem()
    # The final answer is requested in a specific format
    # print(f"\n<<<{result:.4f}>>>")
    # As per instructions, the final result will be at the end.

# Running the function to display the output
solve_kelly_problem()