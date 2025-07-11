import math

def solve_for_F():
    """
    This function calculates the force F required for zero deflection at the beam tip.
    """
    # Step 1: Define geometric parameter 'a'
    a = 12**(1/4)
    a_4 = a**4  # This is simply 12

    # Step 2: Calculate moments of inertia I_ss and I_zz
    # For a square of side H, moment of inertia about centroidal axis is H^4/12
    # Main square (3a x 3a):
    I_main = (3*a)**4 / 12 = 81 * a_4 / 12
    # Cutout square (a x a):
    I_cutout_c = a_4 / 12
    Area_cutout = a**2

    # I_zz: Bending about z-axis. Use parallel axis theorem I = Ic + A*d^2
    # where d is the s-coordinate of the cutout's center.
    d_s = a / 2
    I_zz = I_main - 2 * (I_cutout_c + Area_cutout * d_s**2)
    # I_zz = (81/12)*a_4 - 2*(a_4/12 + a_2 * a_2/4) = (81/12)*a_4 - 2*(a_4/12 + a_4/4) = (81/12)*a_4 - 2*(4*a_4/12)
    # I_zz = (81/12 - 8/12)*a_4 = (73/12)*a_4 = 73

    # I_ss: Bending about s-axis.
    # where d is the z-coordinate of the cutout's center.
    d_z = a
    I_ss = I_main - 2 * (I_cutout_c + Area_cutout * d_z**2)
    # I_ss = (81/12)*a_4 - 2*(a_4/12 + a_2 * a_2) = (81/12)*a_4 - 2*(a_4/12 + a_4) = (81/12)*a_4 - 2*(13*a_4/12)
    # I_ss = (81/12 - 26/12)*a_4 = (55/12)*a_4 = 55


    # Step 3: Calculate L and q0
    L = (30 * I_zz) / 73
    q0 = (9 * I_ss) / 55

    # Step 4: Calculate the force F using the derived formula
    # F = (13 * q0 * L) / 270
    F = (13 * q0 * L) / 270

    # Print the final calculation step-by-step
    print("The final formula to calculate the force F is:")
    print("F = (13 * q0 * L) / 270\n")
    print("Substituting the calculated values for q0 and L:")
    print(f"F = (13 * {q0:.0f} * {L:.0f}) / 270")
    print(f"F = ({13 * q0 * L}) / 270")
    print(f"F = {F}")

solve_for_F()