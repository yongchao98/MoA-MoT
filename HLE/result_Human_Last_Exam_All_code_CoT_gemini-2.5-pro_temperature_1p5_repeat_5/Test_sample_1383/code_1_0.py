import math
from fractions import Fraction

def solve_betting_problem():
    """
    Calculates and explains the optimal and achieved growth rates for a betting scenario.
    """
    p_true_str = ["1/2", "1/4", "1/8", "1/8"]
    q_incorrect_str = ["1/4", "1/2", "1/8", "1/8"]
    
    p_true = [Fraction(s) for s in p_true_str]
    q_incorrect = [Fraction(s) for s in q_incorrect_str]
    odds = [4, 3, 7, 7]

    print("This script solves for the expected logarithmic growth rates in a betting scenario.")
    print("-" * 70)

    # --- Part 1: Calculate W* ---
    # Optimal betting fractions b* are the true probabilities p_true
    b_star = p_true
    r_star = [b * o for b, o in zip(b_star, odds)]

    print("1. Optimal Expected Logarithmic Growth Rate (W*)")
    print("The optimal betting fractions are equal to the true probabilities.")
    print(f"Bets b* = {p_true_str}")
    
    w_star_eq = "W* = "
    terms_star = [f"({p}) * ln({r})" for p, r in zip(p_true, r_star)]
    w_star_eq += " + ".join(terms_star)
    print("The equation for W* is:")
    print(w_star_eq)
    
    # Manually simplified expression for clarity
    print("\nIn simplified form, the expression for W* is:")
    print("W* = (1/2)*ln(2) + (1/4)*ln(3/4) + (1/4)*ln(7/8)")
    print("-" * 70)
    
    # --- Part 2: Calculate W ---
    # The bettor uses fractions b = q_incorrect
    b_w = q_incorrect
    r_w = [b * o for b, o in zip(b_w, odds)]

    print("2. Achieved Growth Rate with Incorrect Probabilities (W)")
    print("The bettor uses fractions equal to their incorrect probability beliefs.")
    print(f"Bets b = {q_incorrect_str}")
    
    w_eq = "W = "
    terms_w = []
    # Note: The expectation is over the true probabilities p_true
    for i in range(len(p_true)):
         terms_w.append(f"({p_true[i]}) * ln({r_w[i]})")
    w_eq += " + ".join(terms_w)
    
    print("The equation for the achieved rate W is:")
    print(w_eq)

    print("\nIn simplified form, the expression for W is:")
    print("W = (1/4)*ln(3/2) + (1/4)*ln(7/8)")
    print("-" * 70)

    # --- Part 3: Calculate Delta W ---
    print("3. Decrease in Growth Rate (ΔW = W* - W)")
    print("The decrease is the difference between the optimal rate and the achieved rate.")
    
    # ΔW is also the KL-Divergence D_KL(P || Q) = Σ p_i * ln(p_i / q_i)
    delta_w_eq = "ΔW = "
    terms_delta = []
    for i in range(len(p_true)):
        if q_incorrect[i] > 0:
            ratio = p_true[i] / q_incorrect[i]
            terms_delta.append(f"({p_true[i]}) * ln({p_true[i]}/{q_incorrect[i]})")
    delta_w_eq += " + ".join(terms_delta)
    
    print("The equation for ΔW is:")
    print(delta_w_eq)

    print("\nThis simplifies to:")
    print("ΔW = (1/2)*ln(2) + (1/4)*ln(1/2)")
    
    print("\nAnd the final result is:")
    print("ΔW = (1/4)*ln(2)")

solve_betting_problem()