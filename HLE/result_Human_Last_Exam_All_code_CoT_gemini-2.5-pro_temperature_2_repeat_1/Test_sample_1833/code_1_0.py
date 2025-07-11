import math

def calculate_tension():
    """
    Calculates the tension in the string for the given physics problem.
    """
    # Step 1: Define the given parameters from the problem
    m = 1.0  # mass of the ring in kg
    M = 1.0  # mass of the object in kg
    g = 9.8  # acceleration due to gravity in m/s^2
    theta_deg = 60.0  # angle of the string in degrees

    # The problem specifies a special case where m = M.
    # The derived formula for tension (T) for this case is:
    # T = (m * g * sin(θ) * (5 + cos(θ)^2)) / (1 + cos(θ)^2)^2
    
    # Step 2: Convert the angle to radians for use in math functions
    theta_rad = math.radians(theta_deg)
    
    # Step 3: Calculate the trigonometric values needed
    sin_theta = math.sin(theta_rad)
    cos_theta_sq = math.cos(theta_rad)**2
    
    # Step 4: Calculate the numerator and denominator of the tension formula
    numerator = m * g * sin_theta * (5 + cos_theta_sq)
    denominator = (1 + cos_theta_sq)**2
    
    # Step 5: Calculate the final tension
    tension = numerator / denominator

    # Step 6: Print the parameters and the final equation with values to show the work
    print("Calculating tension using the formula: T = (m * g * sin(θ) * (5 + cos(θ)²)) / (1 + cos(θ)²)²")
    print(f"\nGiven values:")
    print(f"m = {m} kg")
    print(f"g = {g} m/s²")
    print(f"θ = {theta_deg}°")
    
    print("\nIntermediate values for the equation:")
    print(f"sin({theta_deg}°) = {sin_theta:.4f}")
    print(f"cos({theta_deg}°)^2 = {cos_theta_sq:.4f}")

    print("\nPlugging numbers into the equation:")
    # Using the special case formula (m=M) in the printout for clarity
    print(f"T = ({m:.1f} * {g:.1f} * {sin_theta:.4f} * (5 + {cos_theta_sq:.4f})) / (1 + {cos_theta_sq:.4f})²")
    print(f"T = {numerator:.4f} / {denominator:.4f}")

    # Step 7: Print the final result
    print(f"\nFinal calculated tension: {tension:.4f} N")
    print(f"The tension in the string rounded to two decimal places is {tension:.2f} Newtons.")


calculate_tension()

<<<28.52>>>