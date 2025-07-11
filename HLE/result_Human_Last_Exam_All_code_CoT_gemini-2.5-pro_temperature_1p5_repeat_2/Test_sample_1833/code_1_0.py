import math

def calculate_tension():
    """
    Calculates the tension in the string for the given physics problem.
    """
    # --- Given values ---
    m = 1.0  # kg (mass of the ring)
    M = 1.0  # kg (mass of the object)
    g = 9.8  # m/s^2 (acceleration due to gravity)
    theta_deg = 60.0 # degrees (angle of the string)

    # --- Plan ---
    # The problem is solved using conservation of energy and horizontal momentum,
    # combined with Newton's second law.
    # For the special case where m = M, the derived formula for tension T is:
    # T = g * sin(theta) * (5 + cos(theta)^2) / (1 + cos(theta)^2)^2

    print("Step 1: Define the physical constants and the angle.")
    print(f"Mass of ring m = {m} kg")
    print(f"Mass of object M = {M} kg")
    print(f"Gravity g = {g} m/s^2")
    print(f"Angle theta = {theta_deg} degrees\n")

    # --- Calculations ---
    # Convert angle to radians for use in Python's math functions
    theta_rad = math.radians(theta_deg)

    # Calculate trigonometric values
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)
    cos_sq_theta = cos_theta**2

    # Calculate the tension using the derived formula
    numerator = g * sin_theta * (5 + cos_sq_theta)
    denominator = (1 + cos_sq_theta)**2
    tension = numerator / denominator

    print("Step 2: State the formula for tension (for m=M).")
    print("T = g * sin(theta) * (5 + cos(theta)^2) / (1 + cos(theta)^2)^2\n")

    print("Step 3: Plug in the numbers into the equation.")
    print("First, we evaluate the trigonometric parts:")
    print(f"sin({theta_deg}°) = {sin_theta:.4f}")
    print(f"cos({theta_deg}°) = {cos_theta:.4f}")
    print(f"cos({theta_deg}°)^2 = {cos_sq_theta:.4f}\n")
    
    print("Now, substitute these values into the tension equation:")
    print(f"T = {g} * {sin_theta:.4f} * (5 + {cos_sq_theta:.4f}) / (1 + {cos_sq_theta:.4f})^2")
    print(f"T = {g * sin_theta:.4f} * {5 + cos_sq_theta:.4f} / ({1 + cos_sq_theta:.2f})^2")
    print(f"T = {numerator:.4f} / {denominator:.4f}\n")
    
    print("Step 4: Calculate the final result.")
    print(f"The tension T is {tension:.4f} N.")
    print(f"Rounded to two decimal places, the tension is {tension:.2f} N.")


# Execute the calculation
calculate_tension()