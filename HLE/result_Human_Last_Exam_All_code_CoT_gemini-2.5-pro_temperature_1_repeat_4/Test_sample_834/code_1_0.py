import math

def solve_beam_problem():
    """
    Calculates the force F required for zero deflection at the beam's end.
    """
    print("Step 1: Deriving the formula for F.")
    print("The deflection at x=3L/2 due to F is y_F = 9*F*L^3 / (8*EI).")
    print("The deflection at x=3L/2 due to q(x) is y_q = -13*q0*L^4 / (240*EI).")
    print("Setting y_F + y_q = 0 and solving for F gives: F = (13 * q0 * L) / 270\n")

    print("Step 2: Calculating geometric and load parameters.")
    
    # Given value for a
    a = 12**(1/4)
    a_p4 = 12
    print(f"Given a = 12^(1/4), so a^4 = {a_p4}")

    # Calculate moments of inertia I_ss and I_zz
    # I_ss = I_main - 2 * I_cutout_s = (27*a^4/4) - 2 * (a^4/12 + a^2*a^2) = (81/12 - 2*13/12)*a^4 = 55*a^4/12
    # I_zz = I_main - 2 * I_cutout_z = (27*a^4/4) - 2 * (a^4/12 + a^2*(a/2)^2) = (81/12 - 2*4/12)*a^4 = 73*a^4/12
    I_ss = 55 * a_p4 / 12
    I_zz = 73 * a_p4 / 12
    print(f"Calculated moment of inertia about s-axis, I_ss = {I_ss}")
    print(f"Calculated moment of inertia about z-axis, I_zz = {I_zz}\n")

    # Calculate L and q0 from the given data
    L = 30 * I_zz / 73
    q0 = 9 * I_ss / 55
    print(f"Using the given relations:")
    print(f"L = (30 * I_zz) / 73 = (30 * {I_zz}) / 73 = {L}")
    print(f"q0 = (9 * I_ss) / 55 = (9 * {I_ss}) / 55 = {q0}")
    print("Note: The value of EI is not needed as it cancels out.\n")
    
    print("Step 3: Final calculation of F.")
    # Calculate F using the derived formula and calculated parameters
    F = (13 * q0 * L) / 270
    
    print("Substituting the values into the formula for F:")
    print(f"F = (13 * q0 * L) / 270")
    print(f"F = (13 * {q0} * {L}) / 270")
    print(f"F = {13 * q0 * L} / 270")
    print(f"F = {F}")
    
    return F

if __name__ == "__main__":
    final_force = solve_beam_problem()
    # The final answer is directly returned in the specified format.
    print(f"\n<<< {final_force} >>>")
