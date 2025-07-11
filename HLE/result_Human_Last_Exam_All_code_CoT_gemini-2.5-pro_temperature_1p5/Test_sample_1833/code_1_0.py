import math

def calculate_tension():
    """
    Calculates the tension in the string for the given physics problem.
    """
    # Given values
    m = 1.0  # kg, mass of the ring
    M = 1.0  # kg, mass of the object
    g = 9.8  # m/s^2, acceleration due to gravity
    theta_deg = 60.0  # degrees

    # Convert angle to radians for trigonometric functions
    theta_rad = math.radians(theta_deg)

    # Calculate sin and cos of the angle
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)
    cos_theta_sq = cos_theta**2

    # The derived formula for tension T is:
    # T = (m * M * g * sin(theta) * (3*m + 2*M + M*cos(theta)^2)) / (m + M*cos(theta)^2)^2

    # Calculate the numerator and denominator of the formula
    numerator_val = m * M * g * sin_theta * (3 * m + 2 * M + M * cos_theta_sq)
    denominator_val = (m + M * cos_theta_sq)**2

    # Calculate the tension
    tension = numerator_val / denominator_val

    # Print the equation with the substituted values
    print("The final formula for tension T is:")
    print("T = (m * M * g * sin(theta) * (3*m + 2*M + M*cos(theta)^2)) / (m + M*cos(theta)^2)^2\n")
    print("Plugging in the given values:")
    print(f"m = {m:.1f}, M = {M:.1f}, g = {g:.1f}, theta = {theta_deg:.1f} degrees")
    print(f"sin({theta_deg:.1f}) = {sin_theta:.4f}")
    print(f"cos({theta_deg:.1f})^2 = {cos_theta_sq:.4f}\n")

    print("Calculation:")
    print(f"Numerator = {m:.1f} * {M:.1f} * {g:.1f} * {sin_theta:.4f} * (3*{m:.1f} + 2*{M:.1f} + {M:.1f}*{cos_theta_sq:.4f})")
    print(f"Denominator = ({m:.1f} + {M:.1f} * {cos_theta_sq:.4f})^2")
    print(f"Numerator = {numerator_val:.4f}")
    print(f"Denominator = {denominator_val:.4f}")
    print(f"T = {numerator_val:.4f} / {denominator_val:.4f}\n")

    # Print the final result
    print(f"The calculated tension is {tension:.2f} Newtons.")


calculate_tension()