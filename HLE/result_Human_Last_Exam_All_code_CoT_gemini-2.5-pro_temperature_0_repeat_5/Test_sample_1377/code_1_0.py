import math

def solve_race_growth_rate():
    """
    Calculates the difference between optimal and actual growth rates for a race betting scenario.
    """
    # --- Problem Parameters ---
    # True probabilities for competitors (A, B, C)
    p_A, p_B, p_C = 1/2, 1/4, 1/4
    # Mistaken probabilities
    q_A, q_B, q_C = 1/4, 1/2, 1/4
    # Payout ratios (net odds b:1)
    b_A, b_B, b_C = 4, 3, 3

    # --- 1. Calculate Optimal Growth Rate (W*) ---
    # With true probabilities, only competitor A has a positive edge:
    # Edge_A = p_A * (b_A + 1) - 1 = 1/2 * 5 - 1 = 1.5 > 0
    # Edge_B = p_B * (b_B + 1) - 1 = 1/4 * 4 - 1 = 0
    # Edge_C = p_C * (b_C + 1) - 1 = 1/4 * 4 - 1 = 0
    # Therefore, the optimal strategy is to only bet on A.
    
    # Kelly formula for a single bet: f* = (p*b - (1-p)) / b
    f_A_star = (p_A * b_A - (1 - p_A)) / b_A  # This is 3/8

    # The optimal growth rate W* is the expected log-return.
    # W* = p_A*log(1 + b_A*f_A*) + (1-p_A)*log(1 - f_A*)
    # This simplifies to W* = log(5/4)
    W_star = math.log(5/4)

    # --- 2. Calculate Mistaken Betting Fractions (f_A, f_B) ---
    # With mistaken beliefs, both A and B appear to have a positive edge.
    # The fractions f_A and f_B are found by solving the system of equations:
    # 4*f_A - f_B = q_A*(b_A+1) - 1 = 1/4
    # -f_A + 3*f_B = q_B*(b_B+1) - 1 = 1
    # Solving this system yields f_A = 7/44 and f_B = 17/44.
    f_A_actual = 7/44
    f_B_actual = 17/44

    # --- 3. Calculate Actual Growth Rate (W) ---
    # W is calculated using TRUE probabilities (p) but MISTAKEN fractions (f).
    # W = p_A*log(return if A wins) + p_B*log(return if B wins) + p_C*log(return if C wins)
    
    # Return if A wins: 1 + b_A*f_A - f_B = 1 + 4*(7/44) - 17/44 = 5/4
    # Return if B wins: 1 + b_B*f_B - f_A = 1 + 3*(17/44) - 7/44 = 2
    # Return if C wins: 1 - f_A - f_B = 1 - 7/44 - 17/44 = 20/44 = 5/11
    
    term1_W = p_A * math.log(5/4)
    term2_W = p_B * math.log(2)
    term3_W = p_C * math.log(5/11)
    
    W_actual = term1_W + term2_W + term3_W

    # --- 4. Calculate the Difference (W* - W) ---
    difference = W_star - W_actual

    # Print the final equation with each numerical component.
    print("The final equation is W* - W, calculated as follows:")
    print(f"{W_star:.4f} - ({term1_W:.4f} + {term2_W:.4f} + {term3_W:.4f}) = {difference:.4f}")

solve_race_growth_rate()
<<<0.1354>>>