import math

def solve_race_growth_rates():
    """
    Calculates the difference between optimal (W*) and actual (W) portfolio growth rates
    based on the Kelly Criterion for a three-competitor race.
    """
    # Step 1: Define problem parameters
    # True probabilities of winning for competitors A, B, C
    p_true = [0.5, 0.25, 0.25]
    # Mistakenly believed probabilities
    p_believed = [0.25, 0.5, 0.25]
    # Net odds (b = payout - 1)
    # A has 4:1 payout -> b_A = 3
    # B, C have 3:1 payout -> b_B = 2, b_C = 2
    b = [3, 2, 2]
    competitors = ['A', 'B', 'C']

    print("This problem calculates W* - W, the difference between optimal and actual growth rates.")
    print("-" * 70)

    # Step 2: Calculate the optimal betting fractions (f*) and growth rate (W*)
    print("1. Calculating the Optimal Growth Rate (W*):")
    print(f"   - True Probabilities P(A,B,C): {p_true}")
    print(f"   - Net Odds b(A,B,C): {b}")

    # To find the optimal fractions f*, we first check the expected value of each bet.
    # E_i = p_true_i * b_i - (1 - p_true_i). We only bet if E_i > 0.
    # E_A = 0.5 * 3 - (1 - 0.5) = 1.0.  (Bet on A)
    # E_B = 0.25 * 2 - (1 - 0.25) = -0.25. (Do not bet on B)
    # E_C = 0.25 * 2 - (1 - 0.25) = -0.25. (Do not bet on C)
    # So, f*_B = 0, f*_C = 0. We only need to find f*_A.
    # We maximize G = p_A*log(1+b_A*f_A) + (p_B+p_C)*log(1-f_A).
    # The derivative dG/df_A = 0 gives 3/(1+3f_A) = 1/(1-f_A), which solves to f_A = 1/3.
    f_star = [1/3, 0, 0]
    print(f"   - Optimal betting fractions f*(A,B,C) = [{f_star[0]:.4f}, {f_star[1]:.4f}, {f_star[2]:.4f}]")
    
    # Calculate the wealth factors for each outcome with optimal bets
    w_factor_star_A = 1 - f_star[1] - f_star[2] + b[0] * f_star[0] # A wins
    w_factor_star_B = 1 - f_star[0] - f_star[2] + b[1] * f_star[1] # B wins
    w_factor_star_C = 1 - f_star[0] - f_star[1] + b[2] * f_star[2] # C wins
    
    # Calculate W* using true probabilities and optimal fractions
    w_star = (p_true[0] * math.log(w_factor_star_A) +
              p_true[1] * math.log(w_factor_star_B) +
              p_true[2] * math.log(w_factor_star_C))
    
    print(f"   - The optimal growth rate is W* = P(A)*log({w_factor_star_A:.2f}) + P(B)*log({w_factor_star_B:.2f}) + P(C)*log({w_factor_star_C:.2f})")
    print(f"   - W* = {p_true[0]}*log(2) + {p_true[1]}*log(2/3) + {p_true[2]}*log(2/3) ≈ {w_star:.5f}")
    print("-" * 70)

    # Step 3: Calculate the mistaken betting fractions (f) and actual growth rate (W)
    print("2. Calculating the Actual Growth Rate (W):")
    print(f"   - Believed Probabilities Q(A,B,C): {p_believed}")
    
    # Determine betting fractions 'f' based on mistaken beliefs.
    # E'_A = 0.25 * 3 - (1 - 0.25) = 0. (Do not bet)
    # E'_B = 0.5 * 2 - (1 - 0.5) = 0.5. (Bet)
    # E'_C = 0.25 * 2 - (1 - 0.25) = -0.25. (Do not bet)
    # So, f_A = 0, f_C = 0. We only need to find f_B.
    # dG'/df_B = 0 gives 2/(1+2f_B) = 1/(1-f_B), which solves to f_B = 1/4.
    f_mistaken = [0, 1/4, 0]
    print(f"   - Mistaken betting fractions f(A,B,C) = [{f_mistaken[0]:.4f}, {f_mistaken[1]:.4f}, {f_mistaken[2]:.4f}]")
    
    # Calculate wealth factors with mistaken bets
    w_factor_actual_A = 1 - f_mistaken[1] - f_mistaken[2] + b[0] * f_mistaken[0] # A wins
    w_factor_actual_B = 1 - f_mistaken[0] - f_mistaken[2] + b[1] * f_mistaken[1] # B wins
    w_factor_actual_C = 1 - f_mistaken[0] - f_mistaken[1] + b[2] * f_mistaken[2] # C wins
    
    # Calculate W using true probabilities and mistaken fractions
    w_actual = (p_true[0] * math.log(w_factor_actual_A) +
                p_true[1] * math.log(w_factor_actual_B) +
                p_true[2] * math.log(w_factor_actual_C))

    print(f"   - The actual growth rate is W = P(A)*log({w_factor_actual_A:.2f}) + P(B)*log({w_factor_actual_B:.2f}) + P(C)*log({w_factor_actual_C:.2f})")
    print(f"   - W = {p_true[0]}*log(3/4) + {p_true[1]}*log(3/2) + {p_true[2]}*log(3/4) ≈ {w_actual:.5f}")
    print("-" * 70)

    # Step 4: Calculate the difference
    difference = w_star - w_actual
    print("3. Calculating the Final Difference (W* - W):")
    print(f"   W* - W = (0.5*log(2) + 0.5*log(2/3)) - (0.75*log(3/4) + 0.25*log(3/2))")
    # Expanding the logs to base components log(2) and log(3)
    # W* = log(2) - 0.5*log(3)
    # W = log(3) - 1.75*log(2)
    # W* - W = (log(2) - 0.5*log(3)) - (log(3) - 1.75*log(2))
    print(f"   In terms of elementary logarithms, the equation is:")
    print(f"   W* - W = 2.75 * log(2) - 1.5 * log(3)")
    print(f"   The final result is: {difference:.5f}")


if __name__ == '__main__':
    solve_race_growth_rates()