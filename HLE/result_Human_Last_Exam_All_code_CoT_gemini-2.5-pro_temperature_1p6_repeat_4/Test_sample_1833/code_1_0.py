import math

def calculate_tension():
    """
    Calculates the tension in the string for the given physics problem.
    """
    # Given values from the problem
    m = 1.0  # mass of the ring in kg
    M = 1.0  # mass of the object in kg
    theta_deg = 60.0  # angle in degrees
    g = 9.8  # acceleration due to gravity in m/s^2

    # Convert angle from degrees to radians for use in math functions
    theta_rad = math.radians(theta_deg)

    # Pre-calculate trigonometric values for clarity
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)
    cos_sq_theta = cos_theta**2

    # The formula for tension T is derived from conservation of energy,
    # conservation of horizontal momentum, and Newton's second law.
    # T = (m * M * g * sin(theta) * (3*m + 2*M + M*cos(theta)^2)) / (m + M*cos(theta)^2)^2

    # --- Print the breakdown of the calculation ---
    print("Calculating tension using the formula:")
    print("T = (m * M * g * sin(theta) * (3*m + 2*M + M*cos(theta)^2)) / (m + M*cos(theta)^2)^2")
    print("\nWith the following values:")
    print(f"m = {m} kg")
    print(f"M = {M} kg")
    print(f"g = {g} m/s^2")
    print(f"theta = {theta_deg} degrees\n")

    # Calculate each part of the formula
    numerator_part1 = m * M * g * sin_theta
    numerator_part2 = 3 * m + 2 * M + M * cos_sq_theta
    numerator = numerator_part1 * numerator_part2

    denominator_part = m + M * cos_sq_theta
    denominator = denominator_part**2
    
    print("--- Numerator Calculation ---")
    print(f"m * M * g * sin({theta_deg}°) = {m} * {M} * {g} * {sin_theta:.3f} = {numerator_part1:.3f}")
    print(f"3*m + 2*M + M*cos({theta_deg}°)^2 = 3*{m} + 2*{M} + {M}*({cos_theta:.3f})^2 = {numerator_part2:.3f}")
    print(f"Numerator = {numerator_part1:.3f} * {numerator_part2:.3f} = {numerator:.3f}\n")
    
    print("--- Denominator Calculation ---")
    print(f"m + M*cos({theta_deg}°)^2 = {m} + {M}*({cos_theta:.3f})^2 = {denominator_part:.3f}")
    print(f"Denominator = ({denominator_part:.3f})^2 = {denominator:.3f}\n")
    
    # Calculate the final tension
    tension = numerator / denominator

    print("--- Final Calculation ---")
    print(f"Tension = Numerator / Denominator = {numerator:.3f} / {denominator:.3f}")
    print(f"\nThe calculated tension is {tension:.2f} Newtons.")

# Run the calculation
calculate_tension()