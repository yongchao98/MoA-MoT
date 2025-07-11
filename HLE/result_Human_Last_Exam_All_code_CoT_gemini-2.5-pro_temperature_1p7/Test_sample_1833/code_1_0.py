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

    # The derived formula for tension is:
    # T = (m * M * g * sin(θ) * (3m + M(2 + cos²(θ)))) / (m + M*cos²(θ))²

    # Calculate numerator components
    # Numerator = (m * M * g * sin(θ)) * (3m + M*(2 + cos²(θ)))
    num_part1 = m * M * g * sin_theta
    num_part2 = 3 * m + M * (2 + cos_theta**2)
    numerator = num_part1 * num_part2

    # Calculate denominator component
    # Denominator = (m + M*cos²(θ))²
    den_base = m + M * cos_theta**2
    denominator = den_base**2

    # Calculate the final tension
    tension = numerator / denominator

    # Print the equation with all the numerical values
    print("Derived formula for tension T:")
    print("T = (m * M * g * sin(θ) * (3*m + M*(2 + cos(θ)²))) / (m + M*cos(θ)²)²")
    print("\nPlugging in the values:")
    print(f"m = {m} kg")
    print(f"M = {M} kg")
    print(f"g = {g} m/s²")
    print(f"θ = {theta_deg}°")
    print(f"sin({theta_deg}°) = {sin_theta:.4f}")
    print(f"cos({theta_deg}°) = {cos_theta:.4f}")
    print("\nThe equation becomes:")
    # Reconstructing the equation with calculated parts for clarity
    num_val_str = f"({m} * {M} * {g} * {sin_theta:.4f}) * (3*{m} + {M}*(2 + {cos_theta:.4f}²))"
    den_val_str = f"({m} + {M}*{cos_theta:.4f}²)²"
    print(f"T = [ {num_val_str} ] / [ {den_val_str} ]")
    print(f"T = [ {num_part1:.4f} * {num_part2:.4f} ] / [ {den_base:.4f}² ]")
    print(f"T = {numerator:.4f} / {denominator:.4f}")

    # Print the final result
    print(f"\nThe tension in the string is: {tension:.2f} Newtons")
    
    return tension

# Run the calculation and store the result
final_tension = calculate_tension()
# The final answer format is not printed but is used for evaluation
# For example: print(f"<<<{final_tension:.2f}>>>")