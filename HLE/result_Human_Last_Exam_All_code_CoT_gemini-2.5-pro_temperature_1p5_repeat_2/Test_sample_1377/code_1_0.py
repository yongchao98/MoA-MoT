import math

def solve_kelly_problem():
    """
    Calculates the difference between the optimal growth rate (W*) and the actual
    growth rate (W) based on mistaken probabilities.
    """
    # Step 1: Define variables
    p_true = [0.5, 0.25, 0.25]  # True probabilities for A, B, C
    p_believed = [0.25, 0.5, 0.25] # Believed probabilities for A, B, C
    odds = [3, 2, 2] # Net odds for A, B, C

    # --- Calculation for Optimal Growth Rate (W*) ---
    # Based on true probabilities p_true = (0.5, 0.25, 0.25)
    # The optimal strategy is to bet only on A, as it's the only one with a positive edge.
    # f*_A = p_A - (1-p_A)/o_A = 0.5 - (1-0.5)/3
    f_star_A = p_true[0] - (1 - p_true[0]) / odds[0] # This is 1/3

    # Wealth if A wins (Prob=0.5): 1 + f*_A * o_A = 1 + (1/3)*3 = 2
    # Wealth if A loses (Prob=0.5): 1 - f*_A = 1 - 1/3 = 2/3
    w_star = p_true[0] * math.log(1 + f_star_A * odds[0]) + \
             (p_true[1] + p_true[2]) * math.log(1 - f_star_A)

    # --- Calculation for Actual Growth Rate (W) ---
    # Based on mistaken bets derived from p_believed = (0.25, 0.5, 0.25)
    # The mistaken strategy is to bet only on B, as it's the only one with a positive edge under these beliefs.
    # f'_B = p'_B - (1-p'_B)/o_B = 0.5 - (1-0.5)/2
    f_actual_B = p_believed[1] - (1 - p_believed[1]) / odds[1] # This is 0.25

    # The actual growth rate is calculated using the true probabilities of outcomes.
    # Wealth if A wins (True Prob=0.5): 1 - f'_B = 1 - 0.25 = 0.75
    # Wealth if B wins (True Prob=0.25): 1 + f'_B * o_B = 1 + 0.25*2 = 1.5
    # Wealth if C wins (True Prob=0.25): 1 - f'_B = 1 - 0.25 = 0.75
    w_actual = p_true[0] * math.log(1 - f_actual_B) + \
               p_true[1] * math.log(1 + f_actual_B * odds[1]) + \
               p_true[2] * math.log(1 - f_actual_B)

    # --- Final Calculation ---
    difference = w_star - w_actual

    # The final answer requires outputting the equation components.
    # W* = 0.5*ln(2) + 0.5*ln(2/3) = ln(2) - 0.5*ln(3)
    # W = 0.5*ln(0.75) + 0.25*ln(1.5) + 0.25*ln(0.75) = 0.75*ln(0.75) + 0.25*ln(1.5)
    # 0.75*ln(0.75) = 0.75*(ln(3) - ln(4)) = 0.75*ln(3) - 1.5*ln(2)
    # 0.25*ln(1.5) = 0.25*(ln(3) - ln(2))
    # W = ln(3) - 1.75*ln(2)
    # W* - W = (ln(2) - 0.5*ln(3)) - (ln(3) - 1.75*ln(2))
    #         = 2.75*ln(2) - 1.5*ln(3)
    #         = (11/4)*ln(2) - (3/2)*ln(3)
    ln_2 = math.log(2)
    ln_3 = math.log(3)

    print(f"The optimal growth rate (W*) is calculated as: p_A*ln(1+f*_A*o_A) + (1-p_A)*ln(1-f*_A)")
    print(f"W* = {p_true[0]:.2f}*ln(1 + {f_star_A:.4f}*{odds[0]}) + {(1-p_true[0]):.2f}*ln(1 - {f_star_A:.4f}) = {w_star:.4f}")
    print("\nThe actual growth rate (W) is calculated as: p_A*ln(1-f'_B) + p_B*ln(1+f'_B*o_B) + p_C*ln(1-f'_B)")
    print(f"W = {p_true[0]:.2f}*ln(1-{f_actual_B:.2f}) + {p_true[1]:.2f}*ln(1+{f_actual_B:.2f}*{odds[1]}) + {p_true[2]:.2f}*ln(1-{f_actual_B:.2f}) = {w_actual:.4f}")
    print("\nThe difference W* - W in terms of natural logarithms is:")
    print(f"W* - W = (ln(2) - 0.5*ln(3)) - (ln(3) - 1.75*ln(2))")
    print(f"         = (11/4)*ln(2) - (3/2)*ln(3)")
    print(f"         = (11/4)*{ln_2:.4f} - (3/2)*{ln_3:.4f}")
    print(f"         = {11/4 * ln_2:.4f} - {3/2 * ln_3:.4f}")
    print(f"         = {difference:.4f}")

solve_kelly_problem()