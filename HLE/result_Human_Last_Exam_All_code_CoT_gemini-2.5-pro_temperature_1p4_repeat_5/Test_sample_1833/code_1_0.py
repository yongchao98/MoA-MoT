import math

def calculate_tension():
    """
    Calculates the tension in the string for the given physics problem.
    """
    # 1. Define given constants
    m = 1.0  # kg, mass of the ring
    M = 1.0  # kg, mass of the object
    g = 9.8  # m/s^2, acceleration due to gravity
    theta_deg = 60.0  # degrees

    # 2. Convert angle to radians and get trig values
    theta_rad = math.radians(theta_deg)
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)
    cos2_theta = cos_theta**2

    # 3. Calculate the term L*omega^2 from the energy conservation principle.
    # Formula: L*omega^2 = (2 * g * sin(theta) * (m + M)) / (m + M * cos(theta)^2)
    print("--- Step 1: Calculate the L*omega^2 term ---")
    print("Formula: L*w^2 = (2 * g * sin(theta) * (m + M)) / (m + M * cos(theta)^2)")
    
    num_L_omega_sq = 2 * g * sin_theta * (m + M)
    den_L_omega_sq = m + M * cos2_theta
    L_omega_sq = num_L_omega_sq / den_L_omega_sq

    print(f"L*w^2 = (2 * {g} * {sin_theta:.4f} * ({m} + {M})) / ({m} + {M} * {cos2_theta:.4f})")
    print(f"L*w^2 = {num_L_omega_sq:.4f} / {den_L_omega_sq:.4f} = {L_omega_sq:.4f} m/s^2\n")

    # 4. Calculate the tension T using the force equations.
    # Formula: T = (m * M * (g * sin(theta) + L*omega^2)) / (m + M * cos(theta)^2)
    print("--- Step 2: Calculate the Tension (T) ---")
    print("Formula: T = (m * M * (g * sin(theta) + L*w^2)) / (m + M * cos(theta)^2)")
    
    term_g_sin = g * sin_theta
    num_T = m * M * (term_g_sin + L_omega_sq)
    den_T = m + M * cos2_theta
    T = num_T / den_T

    print(f"T = ({m} * {M} * ({g} * {sin_theta:.4f} + {L_omega_sq:.4f})) / ({m} + {M} * {cos2_theta:.4f})")
    print(f"T = ({m * M} * ({term_g_sin:.4f} + {L_omega_sq:.4f})) / {den_T:.4f}")
    print(f"T = {num_T:.4f} / {den_T:.4f}")
    print(f"T = {T:.4f} N\n")

    # 5. Print the final answer rounded to two decimal places.
    print("--- Final Answer ---")
    print(f"The tension in the string is {T:.2f} Newtons.")
    
    # Return the value for the final answer block
    return T

# Execute the function and capture the result
tension_value = calculate_tension()
# The final answer will be printed inside the function, 
# and also formatted for the platform below.
# print(f'<<<{tension_value:.2f}>>>')

if __name__ == '__main__':
    pass

calculate_tension()
<<<28.52>>>