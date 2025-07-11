from fractions import Fraction

def solve_betting_problem():
    """
    Calculates and explains the achieved growth rate (W) and the decrease
    in growth rate (ΔW) for the given bike racing problem.
    """
    # --- Problem Parameters ---
    # True probabilities
    p = [Fraction(1, 2), Fraction(1, 4), Fraction(1, 8), Fraction(1, 8)]
    # Incorrectly believed probabilities (used for betting fractions)
    q = [Fraction(1, 4), Fraction(1, 2), Fraction(1, 8), Fraction(1, 8)]
    # Odds multipliers (d-for-1)
    d = [4, 3, 7, 7]

    def format_frac(f):
        """Helper function to format fractions for printing."""
        if f.denominator == 1:
            return str(f.numerator)
        return f"{f.numerator}/{f.denominator}"

    # --- Part 1: Calculate the achieved growth rate W ---
    print("--- Achieved Growth Rate (W) ---")
    print("You bet according to your incorrect beliefs (q), so your betting fractions are b_i = q_i.")
    print("The actual growth rate W is calculated with the true probabilities (p_i):\n")
    print("W = p_1*ln(b_1*d_1) + p_2*ln(b_2*d_2) + p_3*ln(b_3*d_3) + p_4*ln(b_4*d_4)")
    
    w_calc_expr = (f"W = ({format_frac(p[0])})*ln({format_frac(q[0])}*{d[0]}) + "
                   f"({format_frac(p[1])})*ln({format_frac(q[1])}*{d[1]}) + "
                   f"({format_frac(p[2])})*ln({format_frac(q[2])}*{d[2]}) + "
                   f"({format_frac(p[3])})*ln({format_frac(q[3])}*{d[3]})")
    print(w_calc_expr)

    w_eval_expr = (f"W = ({format_frac(p[0])})*ln({format_frac(q[0]*d[0])}) + "
                   f"({format_frac(p[1])})*ln({format_frac(q[1]*d[1])}) + "
                   f"({format_frac(p[2])})*ln({format_frac(q[2]*d[2])}) + "
                   f"({format_frac(p[3])})*ln({format_frac(q[3]*d[3]})")
    print(w_eval_expr)
    
    # Since ln(1) = 0, the first term is zero. The last two terms can be combined.
    w_final_expr = f"({format_frac(p[1])})*ln({format_frac(q[1]*d[1])}) + ({format_frac(p[2] + p[3])})*ln({format_frac(q[2]*d[2])})"
    print("\nAfter simplifying (since ln(1)=0):")
    print(f"W = {w_final_expr}\n")

    # --- Part 2: Calculate the decrease in growth rate ΔW ---
    print("--- Decrease in Growth Rate (ΔW) ---")
    print("The decrease is the difference between the optimal rate (W*) and your achieved rate (W).\n")
    
    # Calculate W*
    w_star_final_expr = (f"({format_frac(p[0])})*ln({format_frac(p[0]*d[0])}) + "
                         f"({format_frac(p[1])})*ln({format_frac(p[1]*d[1])}) + "
                         f"({format_frac(p[2] + p[3])})*ln({format_frac(p[2]*d[2])})")
    print("Optimal rate W* (betting b_i = p_i):")
    print(f"W* = {w_star_final_expr}\n")

    # Calculate ΔW = W* - W
    print("The decrease ΔW = W* - W is:")
    print(f"ΔW = [ {w_star_final_expr} ] - [ {w_final_expr} ]")
    
    # Simplify ΔW
    delta_w_unsimplified = (f"({format_frac(p[0])})*ln({format_frac(p[0]*d[0])}) + "
                            f"({format_frac(p[1])})*ln({format_frac(p[1]*d[1])}) - "
                            f"({format_frac(p[1])})*ln({format_frac(q[1]*d[1])})")
    print("\nAfter canceling terms:")
    print(f"ΔW = {delta_w_unsimplified}")

    # Final simplification using log rules
    # This is the KL Divergence D_KL(p||q) = p1*ln(p1/q1) + p2*ln(p2/q2)
    # = (1/2)*ln(2) + (1/4)*ln(1/2) = (1/2)*ln(2) - (1/4)*ln(2) = (1/4)*ln(2)
    delta_w_final = f"({format_frac(p[0] - p[1])})*ln(2)"
    print("\nAfter further simplification:")
    print(f"ΔW = {delta_w_final}")

if __name__ == '__main__':
    solve_betting_problem()