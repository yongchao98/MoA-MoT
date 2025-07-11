import numpy as np
from scipy.integrate import quad

def solve_beam_force():
    """
    Calculates the force F required for zero deflection at the beam's end.
    """

    # Step 1: Define geometric constants and calculate moments of inertia
    a = 12**(1/4)
    a_sq = a**2
    a_p4 = a**4  # This is 12

    print("Step 1: Calculating Moments of Inertia")
    # For a square b x h, I_zc = h*b^3/12, I_sc = b*h^3/12
    # Main square (3a x 3a) centered at (0,0)
    I_zz_main = (3*a) * (3*a)**3 / 12 = 81 * a_p4 / 12
    I_ss_main = (3*a) * (3*a)**3 / 12 = 81 * a_p4 / 12

    # Two cutouts (a x a)
    # Area of one cutout
    A_cutout = a_sq
    # Centroidal moment of inertia of one cutout
    I_zc_cutout = a * a**3 / 12 = a_p4 / 12
    I_sc_cutout = a * a**3 / 12 = a_p4 / 12
    
    # Centers of the cutouts
    # Center 1: (s, z) = (a/2, -a)
    # Center 2: (s, z) = (-a/2, a)
    s_dist_sq_1 = (a/2)**2
    s_dist_sq_2 = (-a/2)**2
    z_dist_sq_1 = (-a)**2
    z_dist_sq_2 = (a)**2

    # Use Parallel Axis Theorem: I = I_c + A*d^2
    # I_zz = I_zz_main - (I_zz_cutout1) - (I_zz_cutout2)
    I_zz_cutout1 = I_zc_cutout + A_cutout * s_dist_sq_1
    I_zz_cutout2 = I_zc_cutout + A_cutout * s_dist_sq_2
    I_zz = I_zz_main - I_zz_cutout1 - I_zz_cutout2

    # I_ss = I_ss_main - (I_ss_cutout1) - (I_ss_cutout2)
    I_ss_cutout1 = I_sc_cutout + A_cutout * z_dist_sq_1
    I_ss_cutout2 = I_sc_cutout + A_cutout * z_dist_sq_2
    I_ss = I_ss_main - I_ss_cutout1 - I_ss_cutout2

    print(f"a^4 = {a_p4}")
    print(f"I_zz = (81/12)a^4 - (a^4/12 + a^2(a/2)^2) - (a^4/12 + a^2(-a/2)^2) = {I_zz}")
    print(f"I_ss = (81/12)a^4 - (a^4/12 + a^2(-a)^2) - (a^4/12 + a^2(a)^2) = {I_ss}")
    print("-" * 20)

    # Step 2: Calculate problem parameters L, q0, and EI
    print("Step 2: Calculating L, q0, and EI")
    L = 30 * I_zz / 73
    q0 = 9 * I_ss / 55
    EI, _ = quad(lambda x: np.sin(x**2), 0, np.pi)

    print(f"L = 30 * I_zz / 73 = 30 * {I_zz:.2f} / 73 = {L:.2f}")
    print(f"q0 = 9 * I_ss / 55 = 9 * {I_ss:.2f} / 55 = {q0:.2f}")
    print(f"EI = integral(sin(x^2)) from 0 to pi = {EI:.4f}")
    print("(Note: EI value is not needed for the final calculation as it cancels out)")
    print("-" * 20)

    # Step 3: Calculate the required force F
    # The condition y_total = 0 gives F = (13 * q0 * L) / 270
    print("Step 3: Calculating Force F")
    F = (13 * q0 * L) / 270

    print("The final equation for F is derived from setting total deflection to zero:")
    print("y_total = y_triangular_load + y_point_force = 0")
    print("13*q0*L^4 / (240*EI) = 9*F*L^3 / (8*EI)")
    print("This simplifies to: F = (13 * q0 * L) / 270")
    print("\nSubstituting the values:")
    print(f"F = (13 * {q0:.2f} * {L:.2f}) / 270")
    print(f"F = {F}")
    print("-" * 20)

    return F

# Execute the function and print the final result
final_F = solve_beam_force()
print("\nFinal Answer:")
print(f"The required force F is: {final_F}")

# The final answer in the required format
print(f"\n<<<{final_F}>>>")
