import math
from fractions import Fraction

def calculate_and_print_growth_rates():
    """
    Calculates and prints the optimal (W*) and actual (W) growth rates,
    and their difference, for the given race betting scenario.
    """
    # Step 1: Define problem parameters
    p_true = [1/2, 1/4, 1/4]
    b = [4, 3, 3]  # Net odds

    # --- W* Calculation (Optimal Strategy) ---
    # Under true probabilities, only betting on A has a positive expectation.
    # The optimal Kelly fraction for a single bet is f* = p - (1-p)/b.
    f_star_A = p_true[0] - (1 - p_true[0]) / b[0]

    # W* is the expected log of wealth with the optimal strategy.
    w_star_outcome_A = 1 + f_star_A * b[0]
    w_star_outcome_loss = 1 - f_star_A

    W_star = (p_true[0] * math.log(w_star_outcome_A) +
              (p_true[1] + p_true[2]) * math.log(w_star_outcome_loss))

    print("--- Optimal Growth Rate (W*) Calculation ---")
    pA_f, pB_f, pC_f = Fraction(p_true[0]), Fraction(p_true[1]), Fraction(p_true[2])
    fA_star_f = Fraction(f_star_A)
    bA_f = Fraction(b[0])
    
    print(f"The optimal strategy is to bet only on A with a fraction of the bankroll f*_A = {fA_star_f}.")
    print("The equation for W* is:")
    print(f"W* = p_A*log(1 + f*_A*b_A) + (p_B + p_C)*log(1 - f*_A)")
    print(f"W* = {pA_f} * log(1 + {fA_star_f} * {bA_f}) + ({pB_f} + {pC_f}) * log(1 - {fA_star_f})")
    print(f"W* = {pA_f} * log({Fraction(w_star_outcome_A)}) + {pB_f + pC_f} * log({Fraction(w_star_outcome_loss)})")
    # W* simplifies to log(5/4)
    print(f"W* = log({Fraction(w_star_outcome_A)**pA_f * Fraction(w_star_outcome_loss)**(pB_f+pC_f)}) = log(5/4) ≈ {W_star:.5f}\n")

    # --- W Calculation (Mistaken Strategy) ---
    # The mistaken fractions are f_A = 1/10, f_B = 3/10, derived from mistaken probabilities.
    f_mistaken = [1/10, 3/10, 0]
    F_mistaken = sum(f_mistaken)

    # W is the expected log of wealth using mistaken fractions with true probabilities.
    w_outcome_A = 1 - F_mistaken + f_mistaken[0] * b[0]
    w_outcome_B = 1 - F_mistaken + f_mistaken[1] * b[1]
    w_outcome_C = 1 - F_mistaken  # No bet on C, so the stake F is lost.

    W = (p_true[0] * math.log(w_outcome_A) +
         p_true[1] * math.log(w_outcome_B) +
         p_true[2] * math.log(w_outcome_C))

    print("--- Actual Growth Rate (W) from Mistaken Strategy ---")
    fA_m_f, fB_m_f = Fraction(f_mistaken[0]), Fraction(f_mistaken[1])
    F_m_f = Fraction(F_mistaken)
    bB_f = Fraction(b[1])
    print(f"The mistaken strategy is to bet f_A = {fA_m_f} and f_B = {fB_m_f}.")
    print("The equation for W is:")
    print(f"W = p_A*log(1-F+f_A*b_A) + p_B*log(1-F+f_B*b_B) + p_C*log(1-F)")
    print(f"W = {pA_f}*log(1-{F_m_f}+{fA_m_f}*{bA_f}) + {pB_f}*log(1-{F_m_f}+{fB_m_f}*{bB_f}) + {pC_f}*log(1-{F_m_f})")
    print(f"W = {pA_f}*log({Fraction(w_outcome_A)}) + {pB_f}*log({Fraction(w_outcome_B)}) + {pC_f}*log({Fraction(w_outcome_C)})")
    # W simplifies to (1/4)*log(9/10) because log(1)=0
    print(f"W = (1/4) * log(9/10) ≈ {W:.5f}\n")
    
    # --- Final Difference Calculation ---
    difference = W_star - W
    print("--- Final Difference (W* - W) ---")
    print("The equation for the difference is:")
    print(f"W* - W = (log(5/4)) - (1/4 * log(9/10))")
    print(f"W* - W = {W_star:.5f} - ({W:.5f})")
    print(f"W* - W = {difference:.5f}")

calculate_and_print_growth_rates()
<<<0.24948>>>