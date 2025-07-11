import math

def calculate_tension():
    """
    Calculates the tension in the string for the given physics problem.
    
    The solution uses the principles of conservation of energy and conservation of
    horizontal momentum to find the angular velocity, and then uses Newton's
    second law in the radial direction to find the tension.
    """
    # Given values
    m = 1.0  # mass of the ring in kg
    M = 1.0  # mass of the object in kg
    g = 9.8  # acceleration due to gravity in m/s^2
    theta_deg = 60.0  # angle in degrees

    # Convert angle to radians for trigonometric functions
    theta_rad = math.radians(theta_deg)
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)

    # Step 1: Find the term L*omega^2 using conservation laws.
    # The formula derived from conservation of energy and horizontal momentum is:
    # L*omega^2 = (2 * g * sin(theta) * (m + M)) / (m + M * cos(theta)^2)
    L_omega_sq_numerator = 2 * g * sin_theta * (m + M)
    L_omega_sq_denominator = m + M * cos_theta**2
    L_omega_sq = L_omega_sq_numerator / L_omega_sq_denominator

    # Step 2: Find the tension T using Newton's second law.
    # The formula for tension, derived from analyzing forces in the radial direction, is:
    # T = (M * g * sin(theta) + M * L*omega^2) / (1 + (M/m) * cos(theta)^2)
    T_numerator = M * g * sin_theta + M * L_omega_sq
    T_denominator = 1 + (M / m) * cos_theta**2
    T = T_numerator / T_denominator

    # Print the detailed calculation steps
    print("--- Calculation Steps ---")
    print("\nStep 1: Calculate the term L*omega^2 (L is string length, omega is angular velocity)")
    print(f"The formula is: L*omega^2 = (2 * g * sin(theta) * (m + M)) / (m + M * cos(theta)^2)")
    print(f"L*omega^2 = (2 * {g} * sin({theta_deg}) * ({m} + {M})) / ({m} + {M} * cos({theta_deg})^2)")
    print(f"L*omega^2 = (2 * {g} * {sin_theta:.4f} * {m + M}) / ({m} + {M} * {cos_theta**2:.4f})")
    print(f"L*omega^2 = {L_omega_sq_numerator:.4f} / {L_omega_sq_denominator:.4f} = {L_omega_sq:.4f} m^2/s^2")

    print("\nStep 2: Calculate the Tension (T)")
    print(f"The formula is: T = (M * g * sin(theta) + M * L*omega^2) / (1 + (M/m) * cos(theta)^2)")
    print(f"T = ({M} * {g} * sin({theta_deg}) + {M} * {L_omega_sq:.4f}) / (1 + ({M}/{m}) * cos({theta_deg})^2)")
    print(f"T = ({M * g * sin_theta:.4f} + {M * L_omega_sq:.4f}) / (1 + {M/m} * {cos_theta**2:.4f})")
    print(f"T = {T_numerator:.4f} / {T_denominator:.4f} = {T:.4f} N")

    # Print the final answer rounded to two decimal places
    print("\n--- Final Answer ---")
    print(f"The tension in the string when theta is {theta_deg} degrees is {T:.2f} Newtons.")
    
    return T

if __name__ == '__main__':
    tension = calculate_tension()
    # The final answer is returned in the required format.
    # print(f"\n<<<{tension:.2f}>>>") # This line is for the final output format.
    final_answer = f"<<<{tension:.2f}>>>"
    # This print is for display, the actual answer is captured in the variable.
    # The final response will only contain the code block and the answer tag.

# Execute the function to generate the output
calculate_tension()