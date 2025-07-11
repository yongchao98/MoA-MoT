import math

def calculate_tension():
    """
    Calculates the tension in the string for the given physics problem.

    The problem involves a ring of mass 'm' sliding on a horizontal rod, with a mass 'M'
    attached by a string of length 'L'. Using conservation of energy, conservation of
    horizontal momentum, and Newton's second law, a formula for the tension (T) can
    be derived as a function of the angle (theta) the string makes with the horizontal.

    For the special case where m = M, the formula simplifies to:
    T = M * g * sin(theta) * (5 + cos(theta)^2) / (1 + cos(theta)^2)^2
    """

    # Define the given constants for the special case
    m = 1.0  # Mass of the ring in kg
    M = 1.0  # Mass of the object in kg
    g = 9.8  # Gravitational acceleration in m/s^2
    theta_deg = 60.0 # Angle in degrees

    # Convert angle to radians for use in math functions
    theta_rad = math.radians(theta_deg)

    # Calculate sine and cosine of the angle
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)
    cos_theta_sq = cos_theta**2

    # Calculate the numerator of the tension formula
    numerator_val = M * g * sin_theta * (5 + cos_theta_sq)

    # Calculate the denominator of the tension formula
    denominator_val = (1 + cos_theta_sq)**2

    # Calculate the final tension
    tension = numerator_val / denominator_val

    # Print the inputs, the formula, and the step-by-step calculation
    print("Problem Parameters:")
    print(f"Mass of ring (m)      = {m} kg")
    print(f"Mass of object (M)    = {M} kg")
    print(f"Angle (theta)         = {theta_deg} degrees")
    print(f"Gravity (g)           = {g} m/s^2\n")

    print("Formula for Tension (T) when m=M:")
    print("T = M * g * sin(theta) * (5 + cos(theta)^2) / (1 + cos(theta)^2)^2\n")

    print("Calculation Breakdown:")
    print(f"sin({theta_deg}°)                  = {sin_theta:.4f}")
    print(f"cos({theta_deg}°)                  = {cos_theta:.4f}")
    print(f"cos({theta_deg}°)^2                 = {cos_theta_sq:.4f}\n")
    
    print("Plugging the numbers into the formula:")
    # Using the variables to show the numbers being used in the final equation
    print(f"Numerator   = {M} * {g} * {sin_theta:.4f} * (5 + {cos_theta_sq:.4f}) = {numerator_val:.4f}")
    print(f"Denominator = (1 + {cos_theta_sq:.4f})^2 = {denominator_val:.4f}\n")

    print(f"T = {numerator_val:.4f} / {denominator_val:.4f}")
    print(f"T = {tension:.4f} N\n")
    
    print("Final Answer:")
    print(f"The tension in the string is {tension:.2f} Newtons.")
    
    # Return the final numerical answer for the platform
    return f"{tension:.2f}"

# Run the calculation and store the final answer
final_answer = calculate_tension()

# The final answer is then wrapped in <<<>>>
# print(f"<<<{final_answer}>>>")