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

    # Convert angle to radians for math functions
    theta_rad = math.radians(theta_deg)

    # Calculate sin and cos of the angle
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)

    # The derived formula for tension T is:
    # T = (m * M * g * sin(theta) * (3m + 2M + M*cos(theta)^2)) / (m + M*cos(theta)^2)^2

    # Break down the calculation for clarity
    numerator_val = m * M * g * sin_theta * (3 * m + 2 * M + M * cos_theta**2)
    denominator_val = (m + M * cos_theta**2)**2

    # Calculate the final tension
    tension = numerator_val / denominator_val
    
    # Print the equation with substituted values for clarity
    print("Step 1: The formula for tension T is:")
    print("T = (m * M * g * sin(θ) * (3*m + 2*M + M*cos(θ)²)) / (m + M*cos(θ)²))²\n")
    
    print(f"Step 2: Substitute the given values (m={m}, M={M}, g={g}, θ={theta_deg}°):")
    # Using formatted strings to show the calculation
    # We use intermediate rounded values for display purposes only
    print(f"sin({theta_deg}°) = {sin_theta:.3f}")
    print(f"cos({theta_deg}°)² = {cos_theta**2:.2f}\n")
    
    num_part1 = f"{m} * {M} * {g} * {sin_theta:.3f}"
    num_part2 = f"(3*{m} + 2*{M} + {M}*{cos_theta**2:.2f})"
    den_part = f"({m} + {M}*{cos_theta**2:.2f})²"
    print(f"T = ({num_part1} * {num_part2}) / {den_part}")
    
    num_calc_val = m * M * g * sin_theta * (3 * m + 2 * M + M * cos_theta**2)
    den_calc_val = (m + M * cos_theta**2)**2
    print(f"T = {num_calc_val:.2f} / {den_calc_val:.4f}\n")

    # Print the final result
    print("Step 3: The calculated tension is:")
    print(f"T = {tension:.2f} Newtons")

# Execute the function
calculate_tension()

<<<28.52>>>