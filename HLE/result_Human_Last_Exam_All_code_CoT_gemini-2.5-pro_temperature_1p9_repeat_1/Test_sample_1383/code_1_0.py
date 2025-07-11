from fractions import Fraction

def solve_betting_problem():
    """
    Solves the bike racing betting problem by calculating optimal and actual
    growth rates and their difference.
    """

    # --- Problem Setup ---
    # True probabilities
    p = [Fraction(1, 2), Fraction(1, 4), Fraction(1, 8), Fraction(1, 8)]
    # Incorrectly believed probabilities
    q = [Fraction(1, 4), Fraction(1, 2), Fraction(1, 8), Fraction(1, 8)]
    # Decimal odds (o-for-1)
    o = [4, 3, 7, 7]
    # Net odds (b = o - 1)
    b = [o_i - 1 for o_i in o]

    print("--- 1. Optimal Growth Rate W* (using true probabilities p) ---")
    f_star = [Fraction(0)] * 4
    for i in range(4):
        edge = p[i] * o[i]
        if edge > 1:
            f_star[i] = (edge - 1) / b[i]
            print(f"Bike {i+1}: p_i * o_i = {p[i]}*{o[i]} = {edge} > 1. This is a favorable bet.")
            print(f"   Optimal fraction f_{i+1}* = ({edge} - 1) / {b[i]} = {f_star[i]}")
        else:
            print(f"Bike {i+1}: p_i * o_i = {p[i]}*{o[i]} = {edge} <= 1. No bet.")
    
    # Symbolic derivation for W*
    # W* = p_1*log(1 - f_1* + f_1*o_1) + (p_2+p_3+p_4)*log(1 - f_1*)
    # W* = (1/2)*log(1 - 1/3 + (1/3)*4) + (1/2)*log(1 - 1/3)
    # W* = (1/2)*log(2) + (1/2)*log(2/3) = (1/2)*log(4/3)
    # W* = log(2) - (1/2)*log(3)
    print("\nThe optimal growth rate W* is (1/2)*log(4/3), which simplifies to log(2) - (1/2)*log(3).\n")
    
    
    print("--- 2. Achieved Growth Rate W (using incorrect beliefs q to bet) ---")
    f = [Fraction(0)] * 4
    for i in range(4):
        edge = q[i] * o[i]
        if edge > 1:
            f[i] = (edge - 1) / b[i]
            print(f"Bike {i+1}: q_i * o_i = {q[i]}*{o[i]} = {edge} > 1. Bettor places a bet.")
            print(f"   Chosen fraction f_{i+1} = ({edge} - 1) / {b[i]} = {f[i]}")
        else:
             print(f"Bike {i+1}: q_i * o_i = {q[i]}*{o[i]} = {edge} <= 1. No bet.")

    # Calculate actual growth rate W using fractions f and true probabilities p
    # W = p_2*log(1 - f_2 + f_2*o_2) + (p_1+p_3+p_4)*log(1 - f_2)
    # W = (1/4)*log(1 - 1/4 + (1/4)*3) + (3/4)*log(1 - 1/4)
    # W = (1/4)*log(6/4) + (3/4)*log(3/4) = (1/4)*log(3/2) + (3/4)*log(3/4)
    # W = (1/4)*(log(3)-log(2)) + (3/4)*(log(3)-log(4))
    # W = (1/4)*log(3)-(1/4)*log(2) + (3/4)*log(3)-(3/4)*2*log(2)
    # W = log(3) - (1/4)*log(2) - (6/4)*log(2) = log(3) - (7/4)*log(2)
    print("\nThe achieved growth rate W, using fractions from incorrect beliefs but with true probabilities, is:")
    print("W = (3/4)*log(3/4) + (1/4)*log(3/2), which simplifies to log(3) - (7/4)*log(2).\n")

    print("--- 3. Decrease in Growth Rate ΔW ---")
    # ΔW = W* - W = (log(2) - (1/2)*log(3)) - (log(3) - (7/4)*log(2))
    # ΔW = log(2) + (7/4)*log(2) - (1/2)*log(3) - log(3)
    # ΔW = (11/4)*log(2) - (3/2)*log(3)
    print("The decrease ΔW = W* - W is calculated as:")
    print("ΔW = (log(2) - (1/2)*log(3)) - (log(3) - (7/4)*log(2))")
    print("This simplifies to (11/4)*log(2) - (3/2)*log(3).")

    # --- Final Answer ---
    w_achieved = "W = log(3) - (7/4)*log(2)"
    delta_w = "ΔW = (11/4)*log(2) - (3/2)*log(3)"

    final_answer = f"""The doubling rate you will achieve is:
{w_achieved}

The amount your doubling rate has decreased is:
{delta_w}"""
    print("\n" + "="*40)
    print("Final Answer:")
    print("="*40)
    print(final_answer)
    print(f"<<<{final_answer}>>>")


if __name__ == "__main__":
    solve_betting_problem()