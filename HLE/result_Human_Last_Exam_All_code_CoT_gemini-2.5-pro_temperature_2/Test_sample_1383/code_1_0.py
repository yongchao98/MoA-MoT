from fractions import Fraction

def solve_betting_problem():
    """
    Calculates the growth rate from betting with incorrect probabilities and
    the decrease from the optimal rate.
    """

    # True probabilities
    p = [Fraction(1, 2), Fraction(1, 4), Fraction(1, 8), Fraction(1, 8)]

    # Incorrectly believed probabilities (used for betting fractions)
    q = [Fraction(1, 4), Fraction(1, 2), Fraction(1, 8), Fraction(1, 8)]

    # Decimal odds (d-for-1 means odds are d)
    d = [4, 3, 7, 7]

    print("--- Part 1: Doubling Rate W from Incorrect Beliefs ---")
    print("The achieved growth rate W is calculated using true probabilities p and betting fractions q:")
    print("W = p1*log(q1*d1) + p2*log(q2*d2) + p3*log(q3*d3) + p4*log(q4*d4)\n")

    # Substitute values into the equation for W
    print("Substituting the given values:")
    term1 = f"({p[0]})*log(({q[0]})*({d[0]}))"
    term2 = f"({p[1]})*log(({q[1]})*({d[1]}))"
    term3 = f"({p[2]})*log(({q[2]})*({d[2]}))"
    term4 = f"({p[3]})*log(({q[3]})*({d[3]}))"
    print(f"W = {term1} + {term2} + {term3} + {term4}\n")

    # Simplify the terms inside the logs
    print("Simplifying the terms inside the logarithms:")
    val1 = q[0] * d[0]
    val2 = q[1] * d[1]
    val3 = q[2] * d[2]
    val4 = q[3] * d[3]
    print(f"W = ({p[0]})*log({val1}) + ({p[1]})*log({val2}) + ({p[2]})*log({val3}) + ({p[3]})*log({val4})\n")

    # Further simplification
    # log(1) is 0, so the first term is 0. The last two terms can be combined.
    print("Final expression for W (since log(1)=0):")
    final_w_expr = f"W = ({p[1]})*log({val2}) + ({p[2] + p[3]})*log({val3})"
    print(final_w_expr)
    print("-" * 50)


    print("\n--- Part 2: Decrease in Doubling Rate ΔW ---")
    print("The decrease ΔW is the KL-divergence D_KL(p || q) = Σ p_i * log(p_i / q_i):")
    print("ΔW = p1*log(p1/q1) + p2*log(p2/q2) + p3*log(p3/q3) + p4*log(p4/q4)\n")

    # Substitute values into the equation for ΔW
    print("Substituting the given values:")
    ratio1 = p[0]/q[0]
    ratio2 = p[1]/q[1]
    ratio3 = p[2]/q[2]
    ratio4 = p[3]/q[3]
    term1_dw = f"({p[0]})*log({p[0]}/{q[0]})"
    term2_dw = f"({p[1]})*log({p[1]}/{q[1]})"
    term3_dw = f"({p[2]})*log({p[2]}/{q[2]})"
    term4_dw = f"({p[3]})*log({p[3]}/{q[3]})"
    print(f"ΔW = {term1_dw} + {term2_dw} + {term3_dw} + {term4_dw}\n")

    # Simplify the terms inside the logs
    print("Simplifying the ratios inside the logarithms:")
    print(f"ΔW = ({p[0]})*log({ratio1}) + ({p[1]})*log({ratio2}) + ({p[2]})*log({ratio3}) + ({p[3]})*log({ratio4})\n")

    # Final simplification using log rules
    # log(1/2) = -log(2) and log(1) = 0
    print("Final expression for ΔW (since log(1)=0 and log(1/2) = -log(2)):")
    print(f"ΔW = ({p[0]})*log(2) + ({p[1]})*(-log(2))")
    print(f"ΔW = (1/2)*log(2) - (1/4)*log(2)")
    final_dw_val = p[0] - p[1]
    print(f"ΔW = ({final_dw_val})*log(2)")
    print("-" * 50)
    
    # Required final answer format
    final_answer = f"W = (1/4)*log(3/2) + (1/4)*log(7/8), ΔW = (1/4)*log(2)"
    print(f"\nFinal Answer: {final_answer}")

if __name__ == "__main__":
    solve_betting_problem()
    final_answer = f"W = (1/4)*log(3/2) + (1/4)*log(7/8), ΔW = (1/4)*log(2)"
    # The user asked me to wrap the answer with "<<<>>>"
    print(f"<<<{final_answer}>>>")