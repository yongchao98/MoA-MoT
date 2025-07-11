import math

def calculate_tension():
    """
    Calculates the tension in the string for the given physics problem.
    """
    # Define the constants from the problem statement
    m = 1.0  # mass of the ring in kg
    M = 1.0  # mass of the object in kg
    g = 9.8  # acceleration due to gravity in m/s^2
    theta_deg = 60.0  # angle in degrees

    # Convert the angle to radians for use in trigonometric functions
    theta_rad = math.radians(theta_deg)

    # Pre-calculate trigonometric values
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)

    # The derived formula for tension (T) is:
    # T = (m * M * g * sin(θ) * (3*m + 2*M + M*cos(θ)²)) / (m + M*cos(θ)²)²

    # Calculate the numerator of the formula
    numerator = (m * M * g * sin_theta * 
                 (3 * m + 2 * M + M * cos_theta**2))
                 
    # Calculate the denominator of the formula
    denominator = (m + M * cos_theta**2)**2

    # Calculate the final tension
    tension = numerator / denominator

    # --- Output the results as requested ---
    
    print("The derived formula for tension T is:")
    print("T = (m*M*g*sin(θ) * (3m + 2M + M*cos²(θ))) / (m + M*cos²(θ))²\n")
    
    print("Plugging in the given values:")
    
    # Print the equation with all the numbers substituted in
    print(f"T = ({m:.1f} * {M:.1f} * {g:.1f} * sin({theta_deg:.0f}°) * (3*{m:.1f} + 2*{M:.1f} + {M:.1f}*cos²({theta_deg:.0f}°))) / ({m:.1f} + {M:.1f}*cos²({theta_deg:.0f}°))²\n")
    
    print(f"Calculating the numerator and denominator:")
    print(f"Numerator = {numerator:.4f}")
    print(f"Denominator = {denominator:.4f}\n")

    print(f"The final calculated tension is: {tension:.2f} N")

calculate_tension()
