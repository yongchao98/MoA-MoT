import math

def calculate_tension():
    """
    Calculates the tension in the string for the given physics problem.
    """
    # 1. Define the given parameters
    M = 1.0  # Mass of the object in kg
    m = 1.0  # Mass of the ring in kg
    g = 9.8  # Acceleration due to gravity in m/s^2
    theta_deg = 60.0 # Angle in degrees

    # 2. Convert angle to radians for trigonometric functions
    theta_rad = math.radians(theta_deg)

    # 3. Calculate the trigonometric components needed for the formula
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)
    cos_sq_theta = cos_theta**2

    # 4. Use the derived formula for tension when m = M:
    # T = M * g * sin(theta) * (5 + cos^2(theta)) / (1 + cos^2(theta))^2
    numerator = M * g * sin_theta * (5 + cos_sq_theta)
    denominator = (1 + cos_sq_theta)**2
    tension = numerator / denominator

    # 5. Print the breakdown of the calculation as requested
    print("Calculating tension T using the formula:")
    print("T = M * g * sin(θ) * (5 + cos²(θ)) / (1 + cos²(θ))²\n")
    print("Plugging in the numbers:")
    print(f"M = {M:.1f} kg")
    print(f"g = {g:.1f} m/s^2")
    print(f"θ = {theta_deg:.1f} degrees")
    print(f"sin({theta_deg}°) = {sin_theta:.4f}")
    print(f"cos²({theta_deg}°) = {cos_sq_theta:.4f}\n")
    
    # Show the final equation with the calculated numbers
    print("The equation becomes:")
    final_equation_str = (
        f"T = {M:.1f} * {g:.1f} * {sin_theta:.4f} * (5 + {cos_sq_theta:.4f}) / "
        f"(1 + {cos_sq_theta:.4f})²"
    )
    print(final_equation_str)
    
    # Show intermediate calculation
    calc_step_2 = f"T = {numerator:.4f} / {denominator:.4f}"
    print(calc_step_2)
    
    # 6. Print the final result rounded to two decimal places
    print(f"\nThe calculated tension is: {tension:.2f} Newtons")

calculate_tension()

# The final numerical answer for the submission system
final_answer = 28.52
print(f"\n<<<{final_answer}>>>")