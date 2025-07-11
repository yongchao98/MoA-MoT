import math

def calculate_tension():
    """
    Calculates the tension in the string for the given physics problem.

    The problem is solved using conservation of energy and momentum to find
    the velocities and accelerations, then Newton's second law to find tension.
    """
    # Given parameters
    m = 1.0  # mass of the ring in kg
    M = 1.0  # mass of the object in kg
    g = 9.8  # acceleration due to gravity in m/s^2
    theta_deg = 60.0
    
    # Convert angle to radians
    theta = math.radians(theta_deg)
    
    # Pre-calculate sine and cosine for convenience
    sin_theta = math.sin(theta)
    cos_theta = math.cos(theta)
    
    # Step 1: Use conservation of energy to find L*w^2 (L=length, w=angular velocity)
    # The derived formula for L*w^2 is: (2 * M * g * sin(theta)) / (M * (m/(m+M)*sin(theta)^2 + cos(theta)^2))
    # This simplifies because L cancels out.
    numerator_Lw2 = 2 * M * g * sin_theta
    denominator_Lw2 = M * ((m / (m + M)) * sin_theta**2 + cos_theta**2)
    L_w_squared = numerator_Lw2 / denominator_Lw2

    # Step 2: Use force analysis in a non-inertial frame to find the ring's acceleration a_m
    # The derived formula for a_m is: (M*cos(theta) / (m + M*cos(theta)^2)) * (L_w_squared + g*sin_theta)
    prefactor_am = (M * cos_theta) / (m + M * cos_theta**2)
    term_am = L_w_squared + g * sin_theta
    a_m = prefactor_am * term_am

    # Step 3: Use the radial force equation to find the tension T
    # T = M*g*sin(theta) - M*a_m*cos(theta) + M*L*w^2
    term1 = M * g * sin_theta
    term2 = -M * a_m * cos_theta
    term3 = M * L_w_squared
    tension = term1 + term2 + term3

    # Print the breakdown of the calculation as requested
    print(f"For the case where m = {m} kg, M = {M} kg, and the angle theta = {theta_deg} degrees:")
    print("\nIntermediate calculated values:")
    print(f"  M*g*sin(theta) = {term1:.2f} N")
    print(f"  -M*a_m*cos(theta) = {term2:.2f} N (contribution from ring's acceleration)")
    print(f"  M*L*w^2 = {term3:.2f} N (centrifugal-like term)")

    print("\nThe tension T is calculated from the sum of forces in the radial direction:")
    print(f"T = M*g*sin(theta) - M*a_m*cos(theta) + M*L*w^2")
    print(f"T = {term1:.2f} N + ({term2:.2f} N) + {term3:.2f} N")
    print(f"T = {tension:.2f} N")


# Run the calculation and print the result
calculate_tension()
# Calculate final value for submission
m, M, g, theta_deg = 1.0, 1.0, 9.8, 60.0
theta = math.radians(theta_deg)
sin_theta = math.sin(theta)
cos_theta = math.cos(theta)
numerator_Lw2 = 2 * M * g * sin_theta
denominator_Lw2 = M * ((m / (m + M)) * sin_theta**2 + cos_theta**2)
L_w_squared = numerator_Lw2 / denominator_Lw2
prefactor_am = (M * cos_theta) / (m + M * cos_theta**2)
term_am = L_w_squared + g * sin_theta
a_m = prefactor_am * term_am
tension = M * g * sin_theta - M * a_m * cos_theta + M * L_w_squared
final_answer = round(tension, 2)
print(f'<<<{final_answer}>>>')
