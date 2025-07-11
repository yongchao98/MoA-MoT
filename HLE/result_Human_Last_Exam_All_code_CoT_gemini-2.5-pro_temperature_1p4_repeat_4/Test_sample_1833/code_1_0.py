import math

def calculate_tension():
    """
    Calculates the tension in the string for the given physics problem.
    """
    # Given values
    m = 1.0  # mass of the ring in kg
    M = 1.0  # mass of the object in kg
    g = 9.8  # acceleration due to gravity in m/s^2
    theta_deg = 60.0  # angle in degrees

    # Convert angle to radians for trigonometric functions
    theta_rad = math.radians(theta_deg)
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)

    # The formula for tension T can be derived from the principles of conservation of energy and momentum.
    # The derived formula is:
    # T = (M * m * g * sin(theta) * (3*m + M*(2 + cos(theta)**2))) / (m + M*cos(theta)**2)**2
    
    # Calculate the components of the formula to show the steps
    cos_theta_sq = cos_theta**2

    # Numerator calculation
    num_part1 = M * m * g * sin_theta
    num_part2 = 3 * m + M * (2 + cos_theta_sq)
    numerator = num_part1 * num_part2

    # Denominator calculation
    den_part1 = m + M * cos_theta_sq
    denominator = den_part1**2
    
    # Final tension calculation
    tension = numerator / denominator

    print("To find the tension, we use a formula derived from conservation laws.")
    print("Formula: T = (M*m*g*sin(theta)*(3*m + M*(2+cos(theta)^2))) / (m + M*cos(theta)^2)^2\n")

    print("Substituting the given values:")
    print(f"Mass of ring (m) = {m} kg")
    print(f"Mass of object (M) = {M} kg")
    print(f"Gravity (g) = {g} m/s^2")
    print(f"Angle (theta) = {theta_deg} degrees\n")

    print("Step-by-step calculation:")
    print("Equation with substituted numbers:")
    print(f"T = ({M}*{m}*{g}*sin({theta_deg}°)*(3*{m} + {M}*(2+cos({theta_deg}°)^2))) / ({m} + {M}*cos({theta_deg}°)^2)^2")
    
    print("\nCalculating each part:")
    print(f"sin(60°) = {sin_theta:.4f}")
    print(f"cos(60°)^2 = {cos_theta_sq:.4f}")

    print(f"\nNumerator = ({M} * {m} * {g} * {sin_theta:.4f}) * (3*{m} + {M}*(2 + {cos_theta_sq:.4f}))")
    print(f"Numerator = ({num_part1:.4f}) * ({num_part2:.4f})")
    print(f"Numerator = {numerator:.4f}")

    print(f"\nDenominator = ({m} + {M} * {cos_theta_sq:.4f})^2")
    print(f"Denominator = ({den_part1:.4f})^2")
    print(f"Denominator = {denominator:.4f}")

    print(f"\nTension T = {numerator:.4f} / {denominator:.4f}")
    
    print(f"\nThe final tension in the string is {tension:.2f} Newtons.")


# Execute the calculation
calculate_tension()