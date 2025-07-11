from fractions import Fraction

def format_frac(f):
    """Helper function to format fractions for printing."""
    if f.denominator == 1:
        return str(f.numerator)
    return f"{f.numerator}/{f.denominator}"

def main():
    """
    Calculates and prints the optimal growth rate (W*), the achieved growth rate (W),
    and the difference (ΔW).
    """
    # True probabilities of each bike winning
    p_true = [Fraction(1, 2), Fraction(1, 4), Fraction(1, 8), Fraction(1, 8)]

    # Incorrectly believed probabilities
    q_belief = [Fraction(1, 4), Fraction(1, 2), Fraction(1, 8), Fraction(1, 8)]

    # Bookmaker odds (d-for-1, so decimal odds are d)
    odds = [Fraction(4), Fraction(3), Fraction(7), Fraction(7)]

    # --- Step 1: Calculate the optimal growth rate W* ---
    # Optimal betting fractions are the true probabilities: b_i = p_i
    b_star = p_true
    
    print("The optimal expected logarithmic growth rate W* is calculated by betting fractions equal to the true probabilities p_i.")
    print("W* = Σ p_i * log(odds_i * p_i)")
    
    # Calculation of W* using properties of logarithms
    # W* = (1/2)log(4*1/2) + (1/4)log(3*1/4) + (1/8)log(7*1/8) + (1/8)log(7*1/8)
    # W* = (1/2)log(2) + (1/4)log(3/4) + (1/4)log(7/8)
    # W* = (1/4) * [2log(2) + log(3/4) + log(7/8)] = (1/4) * log[2^2 * (3/4) * (7/8)] = (1/4) * log(21/8)
    w_star_coeff = Fraction(1, 4)
    w_star_arg = Fraction(21, 8)
    print("The final expression for W* is:")
    print(f"W* = ({format_frac(w_star_coeff)}) * log({format_frac(w_star_arg)})")
    print("-" * 30)

    # --- Step 2: Calculate the achieved growth rate W ---
    # Betting fractions are the incorrect probabilities: b_i = q_i
    b_actual = q_belief
    
    print("The achieved growth rate W is calculated by betting fractions equal to the incorrect probabilities q_i, with the expectation taken over the true probabilities p_i.")
    print("W = Σ p_i * log(odds_i * q_i)")
    
    # Calculation of W
    # W = (1/2)log(4*1/4) + (1/4)log(3*1/2) + (1/8)log(7*1/8) + (1/8)log(7*1/8)
    # W = (1/2)log(1) + (1/4)log(3/2) + (1/4)log(7/8)
    # W = (1/4) * [log(3/2) + log(7/8)] = (1/4) * log[(3/2) * (7/8)] = (1/4) * log(21/16)
    w_coeff = Fraction(1, 4)
    w_arg = Fraction(21, 16)
    print("The final expression for W is:")
    print(f"W = ({format_frac(w_coeff)}) * log({format_frac(w_arg)})")
    print("-" * 30)

    # --- Step 3: Calculate the decrease in growth rate ΔW ---
    # ΔW = W* - W
    print("The decrease in the doubling rate ΔW is the difference between W* and W.")
    print("ΔW = W* - W")

    # Calculation of ΔW
    # ΔW = (1/4)log(21/8) - (1/4)log(21/16)
    # ΔW = (1/4) * log[(21/8) / (21/16)] = (1/4) * log(2)
    delta_w_coeff = Fraction(1, 4)
    delta_w_arg = w_star_arg / w_arg
    print("The final expression for ΔW is:")
    print(f"ΔW = ({format_frac(delta_w_coeff)}) * log({format_frac(delta_w_arg)})")

if __name__ == "__main__":
    main()