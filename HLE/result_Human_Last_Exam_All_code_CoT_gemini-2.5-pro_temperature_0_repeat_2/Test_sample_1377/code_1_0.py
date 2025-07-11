import math
import numpy as np

def solve_growth_rate_difference():
    """
    Calculates the difference between optimal and actual growth rates in a horse race betting scenario.
    """
    # --- Problem Setup ---
    # True probabilities for competitors (A, B, C)
    p_true = np.array([1/2, 1/4, 1/4])
    # Mistaken probabilities
    p_mistaken = np.array([1/4, 1/2, 1/4])
    # Net payout odds (b such that you get b+1 back for a 1 unit bet)
    odds = np.array([4, 3, 3])

    print("--- Step 1: Calculate Optimal Strategy and Growth Rate (W*) ---")
    
    # Check expected profit for each bet with true probabilities
    # E = p*b - (1-p)
    expected_profit_true = p_true * odds - (1 - p_true)
    print(f"Expected profit with true probabilities (A, B, C): {expected_profit_true}")
    
    # For true probabilities, only bet A has positive expectation (1.5). B and C are zero.
    # Thus, the optimal strategy is to only bet on A.
    # The optimal fraction f* for a single bet is p - (1-p)/b
    f_star_A = p_true[0] - (1 - p_true[0]) / odds[0]
    f_star = np.array([f_star_A, 0, 0])
    print(f"Optimal betting fractions f* = (f_A*, f_B*, f_C*) = ({f_star[0]:.4f}, {f_star[1]:.4f}, {f_star[2]:.4f})")

    # Calculate optimal growth rate W*
    # W* = p_A*log(1 + b_A*f_A) + (p_B+p_C)*log(1 - f_A)
    w_star_val = (p_true[0] * math.log(1 + odds[0] * f_star[0]) +
                  p_true[1] * math.log(1 - f_star[0]) +
                  p_true[2] * math.log(1 - f_star[0]))
    
    print(f"Optimal growth rate W* = (1/2)*log(1 + 4*(3/8)) + (1/4)*log(1 - 3/8) + (1/4)*log(1 - 3/8)")
    print(f"W* = log(5/4) = {w_star_val:.5f}\n")

    print("--- Step 2: Calculate Mistaken Betting Fractions (f') ---")
    
    # Check expected profit with mistaken probabilities
    expected_profit_mistaken = p_mistaken * odds - (1 - p_mistaken)
    print(f"Expected profit with mistaken probabilities (A, B, C): {expected_profit_mistaken}")
    
    # With mistaken beliefs, bets A and B have positive expectation, C is zero.
    # We need to solve a system of equations for f'_A and f'_B.
    # The equations derived from maximizing the perceived growth rate are:
    # 37*f_A - 23*f_B = -3
    # 27*f_A +  7*f_B =  7
    A = np.array([[37, -23], [27, 7]])
    B = np.array([-3, 7])
    f_prime_A, f_prime_B = np.linalg.solve(A, B)
    f_prime = np.array([f_prime_A, f_prime_B, 0])
    
    print("Solving the system of equations for the mistaken bettor yields:")
    print(f"Mistaken betting fractions f' = (f'_A, f'_B, f'_C) = ({f_prime[0]:.4f}, {f_prime[1]:.4f}, {f_prime[2]:.4f})")
    print(f"As fractions: f' = (7/44, 17/44, 0)\n")

    print("--- Step 3: Calculate Actual Growth Rate (W) ---")
    
    # Calculate the actual growth rate W using true probabilities and mistaken fractions
    # W = p_A*log(1+b_A*f'_A - f'_B) + p_B*log(1-f'_A+b_B*f'_B) + p_C*log(1-f'_A-f'_B)
    s_A = 1 + odds[0] * f_prime[0] - f_prime[1]
    s_B = 1 - f_prime[0] + odds[1] * f_prime[1]
    s_C = 1 - f_prime[0] - f_prime[1]
    
    w_val = (p_true[0] * math.log(s_A) +
             p_true[1] * math.log(s_B) +
             p_true[2] * math.log(s_C))

    print(f"Actual growth rate W = (1/2)*log(1 + 4*(7/44) - 17/44) + (1/4)*log(1 - 7/44 + 3*(17/44)) + (1/4)*log(1 - 7/44 - 17/44)")
    print(f"W = (1/2)*log(5/4) + (1/4)*log(2) + (1/4)*log(5/11)")
    print(f"W = (1/4)*log(125/88) = {w_val:.5f}\n")

    print("--- Step 4: Calculate the Difference W* - W ---")
    
    difference = w_star_val - w_val
    
    print("The final equation is W* - W:")
    print(f"log(5/4) - (1/4)*log(125/88)")
    print(f"= {w_star_val:.5f} - {w_val:.5f}")
    print(f"= {difference:.5f}")
    
    return difference

if __name__ == '__main__':
    result = solve_growth_rate_difference()
    print(f"\n<<<${result:.5f}>>>")
