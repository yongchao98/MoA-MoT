import math

def calculate_tension():
    """
    Calculates the tension in the string for the given physics problem.
    """
    # Step 1: Define the given physical constants and parameters.
    m = 1.0  # mass of the ring in kg
    M = 1.0  # mass of the object in kg
    g = 9.8  # acceleration due to gravity in m/s^2
    theta_deg = 60.0  # angle in degrees

    # Step 2: Convert the angle to radians for trigonometric functions.
    theta_rad = math.radians(theta_deg)

    # Step 3: Use the derived formula for tension T.
    # T = (M*m*g*sin(theta)*(3*m + 2*M + M*cos^2(theta))) / (m + M*cos^2(theta))^2
    
    # Pre-calculate trigonometric values for clarity
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)
    cos_sq_theta = cos_theta**2

    # Step 4: Calculate the numerator and denominator of the formula.
    numerator = M * m * g * sin_theta * (3 * m + 2 * M + M * cos_sq_theta)
    denominator = (m + M * cos_sq_theta)**2

    # Calculate the final tension
    tension = numerator / denominator

    # Step 5: Print the final equation with all the numerical values substituted in, as requested.
    # This shows how the final answer is derived from the given values.
    print("The tension T is calculated using the formula:")
    print("T = (M * m * g * sin(theta) * (3*m + 2*M + M*cos(theta)^2)) / (m + M*cos(theta)^2)^2")
    print("\nSubstituting the given values:")
    print(f"M = {M} kg, m = {m} kg, g = {g} m/s^2, theta = {theta_deg} degrees")
    # Using f-string to display the equation with numbers plugged in.
    print(f"T = ({M} * {m} * {g} * sin({theta_deg}) * (3*{m} + 2*{M} + {M}*cos({theta_deg})^2)) / ({m} + {M}*cos({theta_deg})^2)^2")

    # Step 6: Print the final calculated tension, rounded to two decimal places.
    print(f"\nThe tension in the string is: {tension:.2f} Newtons")

calculate_tension()
<<<28.52>>>