import math

def solve_and_print():
    """
    This function calculates and prints the symbolic solutions for the doubling rate W
    and the decrease in doubling rate ΔW.
    """

    # Step 1: Define variables symbolically as strings for the final equations.
    # W = p1*log(q1*o1) + p2*log(q2*o2) + p3*log(q3*o3) + p4*log(q4*o4)
    # W = (1/2)*log((1/4)*4) + (1/4)*log((1/2)*3) + (1/8)*log((1/8)*7) + (1/8)*log((1/8)*7)
    # W = (1/2)*log(1) + (1/4)*log(3/2) + (1/8)*log(7/8) + (1/8)*log(7/8)
    # Since log(1) = 0, this simplifies to:
    # W = (1/4)*log(3/2) + (2/8)*log(7/8)
    w_equation_part1_coeff = "1/4"
    w_equation_part1_arg = "3/2"
    w_equation_part2_coeff = "1/4"
    w_equation_part2_arg = "7/8"
    w_final_equation = f"W = ({w_equation_part1_coeff}) * log({w_equation_part1_arg}) + ({w_equation_part2_coeff}) * log({w_equation_part2_arg})"

    # Step 2: Calculate the decrease in doubling rate, ΔW.
    # ΔW = W* - W = D_KL(p || q) = Σ p_i * log(p_i / q_i)
    # ΔW = p1*log(p1/q1) + p2*log(p2/q2) + p3*log(p3/q3) + p4*log(p4/q4)
    # ΔW = (1/2)*log((1/2)/(1/4)) + (1/4)*log((1/4)/(1/2)) + (1/8)*log((1/8)/(1/8)) + (1/8)*log((1/8)/(1/8))
    # ΔW = (1/2)*log(2) + (1/4)*log(1/2) + (1/8)*log(1) + (1/8)*log(1)
    # Since log(1) = 0 and log(1/2) = -log(2), this simplifies to:
    # ΔW = (1/2)*log(2) - (1/4)*log(2)
    delta_w_equation_coeff = "1/4"
    delta_w_equation_arg = "2"
    delta_w_final_equation = f"ΔW = W* - W = ({delta_w_equation_coeff}) * log({delta_w_equation_arg})"

    # Step 3: Print the results in the required format.
    print("The doubling rate W you will achieve is:")
    print(w_final_equation)
    print("\nThe decrease in the doubling rate is:")
    print(delta_w_final_equation)

solve_and_print()