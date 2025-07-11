import sympy
from sympy import init_printing

def display_answer():
    """
    Prints the selected answer choice and the corresponding equations.
    """
    init_printing(use_unicode=True)

    # Define variables
    M0, Rp, R, r, theta = sympy.symbols("M_0 R_p R r theta")
    i_r = sympy.Symbol(r'\hat{i}_r')
    i_theta = sympy.Symbol(r'\hat{i}_\\theta')

    # Expression for H field in Region 1 (0 < r < R_p) from option B
    H1_scalar = M0 * (2*Rp**3 + R**3) / (3*R**3)
    H1_vector = -sympy.cos(theta) * i_r + sympy.sin(theta) * i_theta
    H1 = H1_scalar * H1_vector

    # Expression for H field in Region 2 (R_p < r < R) from option B
    H2_r_scalar = - (2*M0/3) * ((Rp/R)**3 - (Rp/r)**3)
    H2_theta_scalar = (M0/3) * (2*(Rp/R)**3 + (Rp/r)**3)
    H2 = (H2_r_scalar * sympy.cos(theta) * i_r) + (H2_theta_scalar * sympy.sin(theta) * i_theta)

    # Print the answer
    answer_choice = "B"
    print(f"The correct answer is {answer_choice}.\n")
    print("The magnetic field H(r, theta) is given by:\n")
    print("In the region 0 < r < R_p:")
    H1_str_cos = f"H_r = {sympy.pretty(H1_scalar * -sympy.cos(theta), use_unicode=False)}"
    H1_str_sin = f"H_theta = {sympy.pretty(H1_scalar * sympy.sin(theta), use_unicode=False)}"
    print(f"    H = {sympy.pretty(H1, use_unicode=False)}")
    print("Which can be written in components:")
    print(f"    H_r = -M_0 * (2*R_p**3 + R**3)/(3*R**3) * cos(theta)")
    print(f"    H_theta = M_0 * (2*R_p**3 + R**3)/(3*R**3) * sin(theta)\n")


    print("In the region R_p < r < R:")
    H2_str = f"H = {sympy.pretty(H2, use_unicode=False)}"
    print(f"    {H2_str}")
    print("Which can be written in components:")
    print(f"    H_r = {sympy.pretty(H2_r_scalar * sympy.cos(theta), use_unicode=False)}")
    print(f"    H_theta = {sympy.pretty(H2_theta_scalar * sympy.sin(theta), use_unicode=False)}")
    print("\nFinal Answer Selection:")
    print(f"<<<{answer_choice}>>>")

display_answer()