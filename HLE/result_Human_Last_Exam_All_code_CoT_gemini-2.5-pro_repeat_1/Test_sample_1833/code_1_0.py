import math

def calculate_tension():
    """
    Calculates the tension in the string for the given physics problem.
    
    The problem involves a ring of mass m sliding on a horizontal rod,
    with a point mass M attached by a string of length L. The object M
    is released from rest and falls under gravity.

    The calculation proceeds in two main steps:
    1. Use conservation of energy and momentum to find an expression for L*(d(theta)/dt)^2,
       which is related to the centripetal acceleration.
    2. Use Newton's second law in the radial direction to find the tension T.
    """
    
    # Given parameters
    m = 1.0  # mass of the ring in kg
    M = 1.0  # mass of the object in kg
    theta_deg = 60.0  # angle in degrees
    g = 9.8  # acceleration due to gravity in m/s^2

    # Convert angle to radians for trigonometric functions
    theta_rad = math.radians(theta_deg)
    
    # Calculate sin and cos of the angle
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)
    cos_theta_sq = cos_theta**2

    # --- Step 1: Calculate L*omega^2 from energy/momentum conservation ---
    # Formula: L*omega^2 = 2*g*sin(theta)*(m+M) / (m + M*cos^2(theta))
    L_omega_squared_numerator = 2 * g * sin_theta * (m + M)
    L_omega_squared_denominator = m + M * cos_theta_sq
    L_omega_squared = L_omega_squared_numerator / L_omega_squared_denominator

    # --- Step 2: Calculate Tension T ---
    # Formula: T = (m*M / (m + M*cos^2(theta))) * (g*sin(theta) + L*omega^2)
    T_factor = (m * M) / (m + M * cos_theta_sq)
    T_parenthesis = g * sin_theta + L_omega_squared
    tension = T_factor * T_parenthesis

    # --- Output the results ---
    print("This script calculates the tension in the string based on the provided physical setup.")
    print("\n--- Given Values ---")
    print(f"Mass of the ring (m): {m:.2f} kg")
    print(f"Mass of the object (M): {M:.2f} kg")
    print(f"Angle (theta): {theta_deg:.2f} degrees")
    print(f"Gravity (g): {g:.2f} m/s^2")

    print("\n--- Calculation Steps ---")
    print("1. First, we calculate the term L*omega^2 (related to centripetal acceleration):")
    print(f"   L*omega^2 = (2 * {g:.2f} * sin({theta_deg:.2f}°) * ({m:.2f} + {M:.2f})) / ({m:.2f} + {M:.2f} * cos({theta_deg:.2f}°)^2)")
    print(f"   L*omega^2 = {L_omega_squared:.4f} m^2/s^2")
    
    print("\n2. Next, we use this value to calculate the tension (T):")
    print(f"   T = (({m:.2f} * {M:.2f}) / ({m:.2f} + {M:.2f} * cos({theta_deg:.2f}°)^2)) * ({g:.2f} * sin({theta_deg:.2f}°) + {L_omega_squared:.4f})")

    print("\n--- Final Result ---")
    print(f"The tension in the string is {tension:.2f} Newtons.")

calculate_tension()
<<<28.52>>>