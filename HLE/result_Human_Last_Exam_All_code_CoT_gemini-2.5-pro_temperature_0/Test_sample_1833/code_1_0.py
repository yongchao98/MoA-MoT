import math

def calculate_tension():
    """
    Calculates the tension in the string for the given physics problem.
    """
    # Step 1: Define the given values
    m = 1.0  # mass of the ring in kg
    M = 1.0  # mass of the object in kg
    theta_deg = 60.0  # angle in degrees
    g = 9.8  # acceleration due to gravity in m/s^2

    # Step 2: Convert the angle from degrees to radians for trigonometric functions
    theta_rad = math.radians(theta_deg)

    # Step 3: Calculate the sine and cosine of the angle
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)

    # Step 4: Calculate the tension using the derived formula
    # T = (M * m * g * sin(θ) * (3*m + 2*M + M * cos²(θ))) / (m + M * cos²(θ))²
    
    # Calculate the numerator
    numerator = M * m * g * sin_theta * (3 * m + 2 * M + M * cos_theta**2)
    
    # Calculate the denominator
    denominator = (m + M * cos_theta**2)**2
    
    # Calculate the final tension
    tension = numerator / denominator

    # Step 5: Print the final equation with all the numbers substituted
    print("The equation for tension (T) is:")
    print("T = (M * m * g * sin(θ) * (3*m + 2*M + M * cos(θ)²)) / (m + M * cos(θ)²)²")
    print("\nSubstituting the given values:")
    print(f"T = ({M:.1f} * {m:.1f} * {g:.1f} * sin({theta_deg:.1f}°) * (3*{m:.1f} + 2*{M:.1f} + {M:.1f} * cos({theta_deg:.1f}°)²)) / ({m:.1f} + {M:.1f} * cos({theta_deg:.1f}°)²)²")
    print(f"T = ({M:.1f} * {m:.1f} * {g:.1f} * {sin_theta:.4f} * (3*{m:.1f} + 2*{M:.1f} + {M:.1f} * {cos_theta:.4f}²)) / ({m:.1f} + {M:.1f} * {cos_theta:.4f}²)²")
    
    # Step 6: Print the final calculated tension, rounded to two decimal places
    print(f"\nThe calculated tension is {tension:.2f} N.")

# Execute the function
calculate_tension()

<<<28.52>>>