import math

def calculate_tension():
    """
    Calculates the tension in the string for the given physics problem.
    
    The derivation involves using conservation of energy and conservation of
    horizontal momentum to find the velocities, and then Newton's second law
    to find the tension. The final formula for tension T is:
    T = (m * M * g * sin(theta) * (3*m + 2*M + M*cos(theta)^2)) / (m + M*cos(theta)^2)^2
    """
    
    # Given values
    m = 1.0  # mass of the ring in kg
    M = 1.0  # mass of the object in kg
    g = 9.8  # acceleration due to gravity in m/s^2
    theta_deg = 60.0 # angle in degrees
    
    # Convert angle to radians for math functions
    theta_rad = math.radians(theta_deg)
    
    # Calculate trigonometric values
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)
    cos2_theta = cos_theta**2
    
    # Calculate terms in the formula
    numerator_factor1 = m * M * g * sin_theta
    numerator_factor2 = (3 * m) + (2 * M) + (M * cos2_theta)
    denominator = (m + M * cos2_theta)**2
    
    # Calculate the final tension
    tension = (numerator_factor1 * numerator_factor2) / denominator
    
    # Print the breakdown of the calculation as requested
    print("Calculating the tension using the formula:")
    print("T = (m * M * g * sin(theta) * (3*m + 2*M + M*cos(theta)^2)) / (m + M*cos(theta)^2)^2\n")
    print("Plugging in the values:")
    print(f"m = {m} kg")
    print(f"M = {M} kg")
    print(f"g = {g} m/s^2")
    print(f"theta = {theta_deg} degrees")
    print(f"sin(theta) = {sin_theta:.3f}")
    print(f"cos(theta) = {cos_theta:.3f}")
    
    print("\nFinal equation with numbers:")
    # Using f-string to show the plugged-in values in the equation structure
    print(f"T = ({m} * {M} * {g} * {sin_theta:.3f} * (3*{m} + 2*{M} + {M}*{cos2_theta:.3f})) / ({m} + {M}*{cos2_theta:.3f})^2")
    
    # Show intermediate calculation
    num_val = numerator_factor1 * numerator_factor2
    den_val = denominator
    print(f"T = {num_val:.3f} / {den_val:.3f}")

    # Print the final result
    print(f"\nThe tension in the string is: {tension:.2f} N")

calculate_tension()