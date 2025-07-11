import sympy
from sympy import symbols, cos, diff, sqrt, pi, S

def solve_pendulum_period():
    """
    This function symbolically derives the period of the described two-disk pendulum system.
    """
    # Define symbols for mass, radius, and gravity.
    m, R, g = symbols('m R g', positive=True)
    
    # Use a symbol for the angle theta for differentiation
    theta = symbols('theta')

    print("This problem describes a pendulum-like system. We find the period of small oscillations by treating it as a simple harmonic oscillator.")
    print("The period `P` is given by the formula P = 2*pi * sqrt(I_eff / k_eff), where `I_eff` is the effective moment of inertia and `k_eff` is the effective rotational stiffness.\n")

    # --- Step 1: Calculate the effective moment of inertia, I_eff ---
    # For small oscillations, we evaluate the kinetic energy term at equilibrium (theta=0).
    # Let theta_dot be the angular velocity of the rod. I_eff is defined by T = 1/2 * I_eff * theta_dot^2.
    
    # Moment of inertia of a solid disk about its center
    I_disk = S(1)/2 * m * R**2

    # Contributions to I_eff from each component at theta=0:
    # Top Disk (Translational): Conservation of horizontal momentum for the 2m system implies v_1 = -2*R*theta_dot at theta=0.
    # I_eff_1_trans corresponds to m * v_1^2 / theta_dot^2
    I_eff_1_trans = m * (-2*R)**2
    
    # Top Disk (Rotational): Rolling without slipping means omega_1 = v_1/R = -2*theta_dot.
    # I_eff_1_rot corresponds to I_disk * omega_1^2 / theta_dot^2
    I_eff_1_rot = I_disk * (-2)**2
    
    # Bottom Disk (Translational): Its velocity at theta=0 is v_2 = 2*R*theta_dot.
    # I_eff_2_trans corresponds to m * v_2^2 / theta_dot^2
    I_eff_2_trans = m * (2*R)**2
    
    # Bottom Disk (Rotational): Welded to the rod, so omega_2 = theta_dot.
    # I_eff_2_rot corresponds to I_disk * omega_2^2 / theta_dot^2
    I_eff_2_rot = I_disk * (1)**2
    
    # Total effective moment of inertia
    I_eff = sympy.simplify(I_eff_1_trans + I_eff_1_rot + I_eff_2_trans + I_eff_2_rot)

    print(f"1. The effective moment of inertia `I_eff` is calculated from the system's kinetic energy at equilibrium.")
    print(f"   Our calculation yields: I_eff = {I_eff}\n")
    
    # --- Step 2: Calculate the effective spring constant, k_eff ---
    # The potential energy U is due to the lower disk. Let U=0 at the level of the top disk's center.
    # The height of the bottom disk's center is y = -4*R*cos(theta).
    U = m * g * (-4*R*cos(theta))
    
    # k_eff is the second derivative of U with respect to theta, evaluated at equilibrium (theta=0).
    k_eff = sympy.simplify(diff(U, theta, 2).subs(theta, 0))
    
    print(f"2. The effective stiffness `k_eff` is calculated from the system's potential energy.")
    print(f"   Our calculation yields: k_eff = {k_eff}\n")

    # --- Step 3: Calculate the period, P ---
    # P = 2*pi * sqrt(I_eff / k_eff)
    period = 2 * pi * sqrt(I_eff / k_eff)
    
    # Extract numerical coefficients for the explanation
    i_eff_coeff = sympy.simplify(I_eff / (m * R**2))
    k_eff_coeff = sympy.simplify(k_eff / (m * g * R))
    quotient_coeff = sympy.simplify(I_eff / k_eff * g / R)

    print("3. The period of motion is therefore:")
    print(f"   P = 2 * pi * sqrt( I_eff / k_eff )")
    print(f"   P = 2 * pi * sqrt( (({i_eff_coeff}) * m * R**2) / (({k_eff_coeff}) * m * g * R) )")
    print(f"   P = 2 * pi * sqrt( ({quotient_coeff}) * R / g )")
    
    print("\nSo, the final simplified expression for the period is:")
    print(f"   P = {sympy.simplify(period)}")


solve_pendulum_period()