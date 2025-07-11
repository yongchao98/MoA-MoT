import numpy as np
from fractions import Fraction

def solve_growth_rate_difference():
    """
    Calculates the difference between the optimal and actual growth rates
    for a three-competitor race with mistaken probabilities.
    """
    # --- Step 1: Define problem parameters ---
    # True probabilities for competitors (A, B, C)
    p_true = np.array([1/2, 1/4, 1/4])
    # Mistaken probabilities the bettor believes
    p_mistaken = np.array([1/4, 1/2, 1/4])
    # Net odds (b) for a 1 unit bet
    b = np.array([4, 3, 3])
    # Total payout odds (o = b + 1)
    o = b + 1

    # --- Step 2: Calculate Optimal Strategy and Growth Rate (W*) ---
    # We check the expected value for each bet using true probabilities: p_i * o_i
    # E_A = 0.5 * 5 = 2.5 > 1 (Positive expectation, so we bet)
    # E_B = 0.25 * 4 = 1.0 (Zero expectation, so we don't bet)
    # E_C = 0.25 * 4 = 1.0 (Zero expectation, so we don't bet)
    # The problem reduces to a single bet on A.
    # The Kelly fraction for a single bet is f = p - q/b = p - (1-p)/b
    f_A_star = p_true[0] - (1 - p_true[0]) / b[0]
    f_star = np.array([f_A_star, 0, 0])

    # Calculate wealth factors for each outcome under the optimal strategy
    w_star_outcomes = np.array([
        1 + b[0]*f_star[0] - (f_star[1] + f_star[2]), # A wins
        1 + b[1]*f_star[1] - (f_star[0] + f_star[2]), # B wins
        1 + b[2]*f_star[2] - (f_star[0] + f_star[1])  # C wins
    ])

    # Calculate the optimal growth rate W* using true probabilities
    W_star = np.sum(p_true * np.log(w_star_outcomes))

    # --- Step 3: Calculate the Mistaken Strategy (f') ---
    # The bettor calculates fractions based on p_mistaken.
    # Perceived expected value: p'_i * o_i
    # E'_A = 0.25 * 5 = 1.25 > 1 (Bet)
    # E'_B = 0.5 * 4 = 2.0 > 1 (Bet)
    # E'_C = 0.25 * 4 = 1.0 (Don't bet)
    # The bettor will bet on A and B. We solve a system of linear equations
    # derived from maximizing the perceived growth rate G' = p'_A*log(W'_A) + p'_B*log(W'_B) + p'_C*log(W'_C)
    # The system for f'_A and f'_B is:
    # 27*f'_A + 7*f'_B = 7
    # 17*f'_A + 37*f'_B = 17
    A_matrix = np.array([[27, 7], [17, 37]])
    B_vector = np.array([7, 17])
    f_prime_AB = np.linalg.solve(A_matrix, B_vector)
    f_prime = np.array([f_prime_AB[0], f_prime_AB[1], 0])

    # --- Step 4: Calculate the Actual Growth Rate (W) ---
    # The actual growth rate is calculated using the mistaken fractions (f_prime)
    # but with the true probabilities (p_true).
    w_actual_outcomes = np.array([
        1 + b[0]*f_prime[0] - (f_prime[1] + f_prime[2]), # A wins
        1 + b[1]*f_prime[1] - (f_prime[0] + f_prime[2]), # B wins
        1 + b[2]*f_prime[2] - (f_prime[0] + f_prime[1])  # C wins
    ])

    # Calculate the actual growth rate W
    W_actual = np.sum(p_true * np.log(w_actual_outcomes))

    # --- Step 5: Compute and Print the Difference ---
    difference = W_star - W_actual

    print("Calculation of W* - W:")
    print("The difference can be simplified to the expression: (1/4) * log(55/32)")
    print("\nComponent values of this equation:")
    term1 = Fraction(1, 4)
    term2_num = 55
    term2_den = 32
    term2_val = term2_num / term2_den
    log_term2 = np.log(term2_val)

    print(f"First term: {term1.numerator}/{term1.denominator} = {float(term1)}")
    print(f"Argument of log: {term2_num}/{term2_den} = {term2_val}")
    print(f"Natural log of argument: log({term2_val}) = {log_term2:.6f}")
    print(f"\nFinal Equation: W* - W = ({float(term1)}) * ({log_term2:.6f})")
    print(f"Result (W* - W) = {difference:.6f}")
    
    return difference

# Execute the function and print the final answer in the required format
final_answer = solve_growth_rate_difference()
print(f"\n<<<{final_answer:.6f}>>>")
