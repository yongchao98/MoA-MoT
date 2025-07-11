import math

def calculate_tension():
    """
    Calculates the tension in the string for the given physics problem.
    """
    # Given parameters from the problem description
    m = 1.0  # kg, mass of the ring
    M = 1.0  # kg, mass of the object
    g = 9.8  # m/s^2, acceleration due to gravity
    theta_deg = 60.0  # degrees, the angle of the string

    # --- Step 1: Preliminary calculations ---
    # Convert the angle from degrees to radians for use in trigonometric functions
    theta_rad = math.radians(theta_deg)

    # Calculate the sine and cosine of the angle
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)

    # --- Step 2: Calculate L*omega^2 from energy conservation ---
    # The term L*omega^2 is derived from the conservation of energy principle.
    # Formula: L*omega^2 = (2*g*sin(theta)) / [ (m*sin(theta)^2)/(m+M) + cos(theta)^2 ]
    l_omega_sq_numerator = 2 * g * sin_theta
    l_omega_sq_denominator = (m * sin_theta**2) / (m + M) + cos_theta**2
    l_omega_sq = l_omega_sq_numerator / l_omega_sq_denominator

    # --- Step 3: Calculate the tension (T) using Newton's second law ---
    # The tension is found from the force balance on the mass M.
    # Formula: T = (M*m*(L*omega^2 + g*sin(theta))) / (m + M*cos(theta)^2)
    tension_numerator = M * m * (l_omega_sq + g * sin_theta)
    tension_denominator = m + M * cos_theta**2
    tension = tension_numerator / tension_denominator

    # --- Step 4: Output the results with explanation ---
    print("To solve the problem, we first find the term L*omega^2 and then use it to find the tension T.")
    
    print("\nCalculation of L*omega^2:")
    print(f"L*omega^2 = (2 * {g:.1f} * sin({theta_deg:.0f}°)) / (({m:.1f} * sin({theta_deg:.0f}°)^2) / ({m:.1f} + {M:.1f}) + cos({theta_deg:.0f}°)^2)")
    print(f"L*omega^2 = (2 * {g:.1f} * {sin_theta:.2f}) / (({m:.1f} * {sin_theta**2:.2f}) / {m + M:.1f} + {cos_theta**2:.2f})")
    print(f"L*omega^2 = {l_omega_sq_numerator:.2f} / {l_omega_sq_denominator:.2f} = {l_omega_sq:.2f} m/s^2")

    print("\nCalculation of Tension (T):")
    print(f"T = ({M:.1f} * {m:.1f} * (L*omega^2 + {g:.1f} * sin({theta_deg:.0f}°))) / ({m:.1f} + {M:.1f} * cos({theta_deg:.0f}°)^2)")
    print(f"T = ({M:.1f} * {m:.1f} * ({l_omega_sq:.2f} + {g:.1f} * {sin_theta:.2f})) / ({m:.1f} + {M:.1f} * {cos_theta**2:.2f})")
    print(f"T = {tension_numerator:.2f} / {tension_denominator:.2f}")

    print(f"\nThe final tension in the string is {tension:.2f} Newtons.")
    
    # Final answer in the specified format
    print(f"<<<{tension:.2f}>>>")

# Execute the function
calculate_tension()