from fractions import Fraction

def solve_betting_problem():
    """
    This script calculates the achieved logarithmic growth rate (W) and
    the decrease from the optimal rate (Delta W) for a betting scenario
    with incorrect probability beliefs.
    """

    # --- Problem Parameters ---
    # p: true probabilities
    # q: incorrect (believed) probabilities
    # o: decimal odds (k-for-1 means odds are k)
    p = [Fraction(1, 2), Fraction(1, 4), Fraction(1, 8), Fraction(1, 8)]
    q = [Fraction(1, 4), Fraction(1, 2), Fraction(1, 8), Fraction(1, 8)]
    o = [4, 3, 7, 7]
    
    # --- Part 1: Calculate the achieved growth rate W ---
    
    print("1. Calculating the achieved growth rate W")
    print("The betting fractions 'b' are set to the believed probabilities 'q'.")
    print("The achieved growth rate W is calculated with true probabilities 'p':")
    print("W = p_1*ln(q_1*o_1) + p_2*ln(q_2*o_2) + p_3*ln(q_3*o_3) + p_4*ln(q_4*o_4)")
    print(f"W = ({p[0]})*ln({q[0]} * {o[0]}) + ({p[1]})*ln({q[1]} * {o[1]}) + ({p[2]})*ln({q[2]} * {o[2]}) + ({p[3]})*ln({q[3]} * {o[3]})")
    print("Simplifying the terms inside the logs:")
    print(f"W = ({p[0]})*ln({q[0]*o[0]}) + ({p[1]})*ln({q[1]*o[1]}) + ({p[2]})*ln({q[2]*o[2]}) + ({p[3]})*ln({q[3]*o[3]})")
    print("Using ln(a*b) = ln(a) + ln(b) and simplifying further:")
    print("W = (1/2)*ln(1) + (1/4)*(ln(3) - ln(2)) + (1/8)*(ln(7) - ln(8)) + (1/8)*(ln(7) - ln(8))")
    print("Using ln(1)=0 and ln(8)=3*ln(2):")
    print("W = 0 + (1/4)*ln(3) - (1/4)*ln(2) + (1/8)*ln(7) - (3/8)*ln(2) + (1/8)*ln(7) - (3/8)*ln(2)")
    print("Combining coefficients for ln(2), ln(3), and ln(7):")
    
    w_coeff_ln2 = Fraction(-1, 4) - Fraction(3, 8) - Fraction(3, 8)
    w_coeff_ln3 = Fraction(1, 4)
    w_coeff_ln7 = Fraction(1, 8) + Fraction(1, 8)
    
    print(f"W = ({w_coeff_ln3})*ln(3) + ({w_coeff_ln7})*ln(7) {w_coeff_ln2}*ln(2)")
    print("\nThe final expression for W is:")
    print(f"W = {w_coeff_ln3}*ln(3) + {w_coeff_ln7}*ln(7) - {abs(w_coeff_ln2)}*ln(2)")
    
    print("\n" + "="*50 + "\n")

    # --- Part 2: Calculate the decrease in growth rate Delta W ---

    print("2. Calculating the decrease in growth rate Delta W")
    print("Delta W = W* - W, where W* is the optimal rate from betting b=p.")
    print("This simplifies to the KL-divergence D(p||q):")
    print("Delta W = p_1*ln(p_1/q_1) + p_2*ln(p_2/q_2) + p_3*ln(p_3/q_3) + p_4*ln(p_4/q_4)")
    print(f"Delta W = ({p[0]})*ln(({p[0]})/({q[0]})) + ({p[1]})*ln(({p[1]})/({q[1]})) + ({p[2]})*ln(({p[2]})/({q[2]})) + ({p[3]})*ln(({p[3]})/({q[3]}))")
    print("Simplifying the ratios:")
    print(f"Delta W = ({p[0]})*ln({p[0]/q[0]}) + ({p[1]})*ln({p[1]/q[1]}) + ({p[2]})*ln({p[2]/q[2]}) + ({p[3]})*ln({p[3]/q[3]})")
    print("Using ln(1/2) = -ln(2) and ln(1) = 0:")
    print("Delta W = (1/2)*ln(2) + (1/4)*(-ln(2)) + 0 + 0")
    print("Combining coefficients for ln(2):")
    
    delta_w_coeff_ln2 = p[0] - p[1]
    
    print(f"Delta W = ({delta_w_coeff_ln2})*ln(2)")
    print("\nThe final expression for Delta W is:")
    print(f"Delta W = {delta_w_coeff_ln2}*ln(2)")

if __name__ == '__main__':
    solve_betting_problem()