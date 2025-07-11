import math

def calculate_tension():
    """
    Calculates the tension in the string for the given physics problem.
    """
    # --- Given values ---
    m = 1.0  # mass of the ring in kg
    M = 1.0  # mass of the object in kg
    g = 9.8  # gravitational acceleration in m/s^2
    theta_deg = 60.0  # angle in degrees

    # Convert angle to radians for calculations
    theta_rad = math.radians(theta_deg)
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)

    # For the special case m=M, the equations simplify.
    
    # Step 1: Calculate the term L*(d(theta)/dt)^2 from energy conservation.
    # This term represents the square of the relative speed between the masses.
    # L_omega_sq stands for L * omega^2
    denominator_L_omega_sq = 1 - 0.5 * sin_theta**2
    L_omega_sq = (2 * g * sin_theta) / denominator_L_omega_sq

    # Step 2: To find the ring's acceleration a_rx, we first need L*(d^2(theta)/dt^2).
    # This can be found by solving a system of equations derived from the tangential force equation.
    # For m=M, a_rx = (L/2)*(cos(theta)*omega^2 + sin(theta)*theta_ddot)
    # and L*theta_ddot = g*cos(theta) + a_rx*sin(theta)
    # Solving for a_rx:
    a_rx_numerator = (L_omega_sq / 2) * cos_theta
    a_rx_denominator = 1 - (sin_theta**2 / 2)
    a_rx = a_rx_numerator / a_rx_denominator

    # Step 3: Calculate the tension T using the radial force equation.
    # T = M*g*sin(theta) - M*a_rx*cos(theta) + M*L*(omega^2)
    term1_gravity = M * g * sin_theta
    term2_inertial = -M * a_rx * cos_theta
    term3_centripetal = M * L_omega_sq
    
    T = term1_gravity + term2_inertial + term3_centripetal
    
    # Output the final equation with numerical values
    print("The final tension T is the sum of three components:")
    print("1. Component from gravity: M*g*sin(theta)")
    print("2. Component from the ring's acceleration (inertial force): -M*a_rx*cos(theta)")
    print("3. Component from the relative circular motion (centripetal force): M*L*(d(theta)/dt)^2\n")

    print("Final Equation with calculated values:")
    print(f"T = ({M:.1f} * {g:.1f} * sin({theta_deg:.0f})) - ({M:.1f} * {a_rx:.2f} * cos({theta_deg:.0f})) + ({M:.1f} * {L_omega_sq:.2f})")
    print(f"T = {term1_gravity:.2f} N + ({term2_inertial:.2f} N) + {term3_centripetal:.2f} N")
    print(f"T = {T:.2f} N")

calculate_tension()