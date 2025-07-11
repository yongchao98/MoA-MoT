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

    # Convert angle to radians for use in trigonometric functions
    theta_rad = math.radians(theta_deg)

    # Calculate sin(theta) and cos(theta)
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)
    cos_theta_sq = cos_theta**2

    # The derived formula for tension is:
    # T = (m * M * g * sin(theta) * (3*m + 2*M + M*cos(theta)^2)) / (m + M*cos(theta)^2)^2

    # Calculate the numerator and denominator separately for clarity
    numerator_val = (m * M * g * sin_theta) * (3 * m + 2 * M + M * cos_theta_sq)
    denominator_val = (m + M * cos_theta_sq)**2

    # Calculate the final tension
    tension = numerator_val / denominator_val

    # --- Output ---
    print("The derived formula for tension (T) is:")
    print("T = (m * M * g * sin(theta) * (3*m + 2*M + M*cos(theta)^2)) / (m + M*cos(theta)^2)^2")
    print("\nSubstituting the given values:")
    print(f"m = {m} kg")
    print(f"M = {M} kg")
    print(f"g = {g} m/s^2")
    print(f"theta = {theta_deg} degrees")
    print("\nThe equation with each number substituted is:")
    
    # Building the string for the equation with substituted values
    num_part1_str = f"{m} * {M} * {g} * {sin_theta:.4f}"
    num_part2_str = f"(3*{m} + 2*{M} + {M}*({cos_theta:.4f})^2)"
    den_part_str = f"({m} + {M}*({cos_theta:.4f})^2)^2"
    
    print(f"T = ({num_part1_str}) * {num_part2_str} / {den_part_str}")

    print(f"\nCalculated tension: {tension:.2f} Newtons")


calculate_tension()
<<<28.52>>>