from fractions import Fraction

def format_ln_expression(coeff2, coeff3):
    """Formats an expression of the form c2*ln(2) + c3*ln(3) into a readable string."""
    parts = []
    # Format ln(2) term
    if coeff2 != 0:
        if abs(coeff2) == 1:
            term = "ln(2)" if coeff2 > 0 else "-ln(2)"
        else:
            term = f"{coeff2}*ln(2)" if coeff2 > 0 else f"-{abs(coeff2)}*ln(2)"
        parts.append(term)
    
    # Format ln(3) term
    if coeff3 != 0:
        sign = " + " if coeff3 > 0 else " - "
        
        if abs(coeff3) == 1:
            term = "ln(3)"
        else:
            term = f"{abs(coeff3)}*ln(3)"
        
        if not parts: # if ln(2) term was zero
            if coeff3 > 0:
                parts.append(term)
            else:
                parts.append(f"-{term}")
        else:
            parts.append(f"{sign}{term}")

    if not parts:
        return "0"
    
    # Capitalize the first letter if it's not a minus sign
    result = "".join(parts)
    if result.startswith('-'):
        return result
    else:
        # Remove leading " + " if it exists (e.g. from a positive ln(3) with zero ln(2))
        return result.lstrip(" + ")


def solve_and_print():
    """
    Solves the bike race betting problem and prints the results.
    """
    p = [Fraction(1, 2), Fraction(1, 4), Fraction(1, 8), Fraction(1, 8)]
    q = [Fraction(1, 4), Fraction(1, 2), Fraction(1, 8), Fraction(1, 8)]
    o = [4, 3, 7, 7]

    # --- Part 1: Calculate W* ---
    # Find favorable bets for true probabilities p: p_i * o_i > 1
    # p1*o1 = 1/2*4=2 > 1. Bet on Bike 1.
    # p2*o2=1/4*3=3/4<1, p3*o3=1/8*7=7/8<1, p4*o4=1/8*7=7/8<1.
    b1_star = (p[0] * o[0] - 1) / (o[0] - 1)
    
    # Calculate returns for W*
    ret_win1 = 1 - b1_star + b1_star * o[0]
    ret_lose1 = 1 - b1_star
    
    # --- Part 2: Calculate W ---
    # Find favorable bets for incorrect probabilities q: q_i * o_i > 1
    # q1*o1=1/4*4=1, q2*o2=1/2*3=3/2>1. Bet on Bike 2.
    # q3*o3=7/8<1, q4*o4=7/8<1.
    b2_prime = (q[1] * o[1] - 1) / (o[1] - 1)

    # Calculate returns for W, using true probabilities p
    ret_win2 = 1 - b2_prime + b2_prime * o[1]
    ret_lose2 = 1 - b2_prime

    # --- Calculations for printing ---
    # W* = p1*ln(ret_win1) + p2*ln(ret_lose1) + p3*ln(ret_lose1) + p4*ln(ret_lose1)
    # W* = 1/2*ln(2) + 1/4*ln(2/3) + 1/8*ln(2/3) + 1/8*ln(2/3)
    c2_w_star = p[0] * 1 + (p[1]+p[2]+p[3]) * 1 # ln(2) from ln(2) and ln(2/3)
    c3_w_star = p[0] * 0 + (p[1]+p[2]+p[3]) * (-1) # ln(3) from ln(2/3)
    w_star_final = format_ln_expression(c2_w_star, c3_w_star)

    # W = p1*ln(ret_lose2) + p2*ln(ret_win2) + p3*ln(ret_lose2) + p4*ln(ret_lose2)
    # W = 1/2*ln(3/4) + 1/4*ln(3/2) + 1/8*ln(3/4) + 1/8*ln(3/4)
    # Decomposing ln(3/2) = ln(3)-ln(2); ln(3/4) = ln(3)-2ln(2)
    c2_w = p[1]*(-1) + (p[0]+p[2]+p[3])*(-2) # ln(2) coeff
    c3_w = p[1]*(1) + (p[0]+p[2]+p[3])*(1) # ln(3) coeff
    w_final = format_ln_expression(c2_w, c3_w)

    # --- Delta W ---
    delta_w_c2 = c2_w_star - c2_w
    delta_w_c3 = c3_w_star - c3_w
    delta_w_final = format_ln_expression(delta_w_c2, delta_w_c3)

    # --- Print results ---
    print("The doubling rate achieved with the incorrect probabilities, W, is:")
    print(f"W = ({p[0]})*ln({ret_lose2}) + ({p[1]})*ln({ret_win2}) + ({p[2]})*ln({ret_lose2}) + ({p[3]})*ln({ret_lose2})")
    print(f"W = {w_final}\n")
    
    print("The optimal expected doubling rate, W*, is:")
    print(f"W* = ({p[0]})*ln({ret_win1}) + ({p[1]})*ln({ret_lose1}) + ({p[2]})*ln({ret_lose1}) + ({p[3]})*ln({ret_lose1})")
    print(f"W* = {w_star_final}\n")

    print("The decrease in the doubling rate, Delta W = W* - W, is:")
    print(f"Delta W = ({w_star_final}) - ({w_final})")
    print(f"Delta W = {delta_w_final}")


solve_and_print()
<<<W = ln(3) - 7/4*ln(2), Delta W = 11/4*ln(2) - 3/2*ln(3)>>>