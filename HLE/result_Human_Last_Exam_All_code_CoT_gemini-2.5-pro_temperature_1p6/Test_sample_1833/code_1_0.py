import math

def calculate_tension():
    """
    Calculates the tension in the string for the given physics problem.
    """
    # Given values
    m = 1.0  # Mass of the ring in kg
    M = 1.0  # Mass of the object in kg
    g = 9.8  # Acceleration due to gravity in m/s^2
    theta_deg = 60.0  # Angle in degrees

    # Convert angle to radians for trigonometric functions
    theta_rad = math.radians(theta_deg)
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)

    # The derived formula for tension T is:
    # T = (m * M * g * sin(theta) * (3*m + 2*M + M * cos(theta)^2)) / (m + M * cos(theta)^2)^2
    # For the special case m = M, this simplifies to:
    # T = M * g * sin(theta) * (5 + cos(theta)^2) / (1 + cos(theta)^2)^2
    
    # Calculate numerator and denominator of the simplified formula
    numerator_val = M * g * sin_theta * (5 + cos_theta**2)
    denominator_val = (1 + cos_theta**2)**2
    
    # Calculate the tension
    tension = numerator_val / denominator_val

    # Output the steps and the final answer
    print("The formula for tension T when m = M is:")
    print("T = M * g * sin(theta) * (5 + cos(theta)^2) / (1 + cos(theta)^2)^2\n")
    print("Plugging in the values:")
    print(f"M = {M} kg")
    print(f"g = {g} m/s^2")
    print(f"theta = {theta_deg} degrees\n")

    print("The equation with the numbers substituted is:")
    # Show the numbers in the equation
    print(f"T = {M} * {g} * sin({theta_deg}) * (5 + (cos({theta_deg}))^2) / (1 + (cos({theta_deg}))^2)^2")
    print(f"T = {M} * {g} * {sin_theta:.4f} * (5 + {cos_theta**2:.4f}) / (1 + {cos_theta**2:.4f})^2")
    print(f"T = {M*g*sin_theta:.4f} * {5 + cos_theta**2:.4f} / {1 + cos_theta**2:.4f}^2")
    print(f"T = {numerator_val:.4f} / {denominator_val:.4f}\n")
    
    # Final Answer
    print(f"The tension in the string is {tension:.2f} Newtons.")


# Run the calculation
calculate_tension()

print("\n<<<28.52>>>")