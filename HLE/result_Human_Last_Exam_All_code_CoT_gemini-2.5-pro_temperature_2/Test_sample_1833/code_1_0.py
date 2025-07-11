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

    # Since m=M, the general formula for tension T simplifies.
    # The derived formula for T when m=M is:
    # T = g * sin(θ) * (5 + cos²(θ)) / (1 + cos²(θ))²

    # Convert angle to radians for use in math functions
    theta_rad = math.radians(theta_deg)

    # Calculate sine and cosine of the angle
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)
    cos_sq_theta = cos_theta**2

    # Calculate the numerator and denominator of the main fraction
    numerator = 5 + cos_sq_theta
    denominator_base = 1 + cos_sq_theta
    denominator = denominator_base**2
    
    # Calculate the tension
    tension = g * sin_theta * (numerator / denominator)

    print("Step-by-step calculation for the tension T:")
    print("Given formula for the special case m = M = 1kg:")
    print(f"T = g * sin(θ) * (5 + cos²(θ)) / (1 + cos²(θ))²\n")
    print("Plugging in the values:")
    print(f"g = {g} m/s²")
    print(f"θ = {theta_deg}°\n")
    
    # Show the calculation with intermediate numerical values
    print(f"1. Calculate components:")
    print(f"   sin({theta_deg}°) = {sin_theta:.4f}")
    print(f"   cos²({theta_deg}°) = {cos_sq_theta:.4f}\n")
    
    print(f"2. Substitute into the formula:")
    print(f"   T = {g} * {sin_theta:.4f} * (5 + {cos_sq_theta:.4f}) / (1 + {cos_sq_theta:.4f})²")
    print(f"   T = {g} * {sin_theta:.4f} * ({numerator:.4f}) / ({denominator_base:.4f})²")
    print(f"   T = {g} * {sin_theta:.4f} * ({numerator:.4f}) / ({denominator:.4f})")
    print(f"   T = {g} * {sin_theta:.4f} * {numerator / denominator:.4f}\n")
    
    print(f"3. Final Calculation:")
    print(f"   T = {tension:.4f} N\n")

    # Print the final answer rounded to two decimal places
    print("The final tension in the string, rounded to two decimal places, is:")
    print(f"{tension:.2f} N")

calculate_tension()