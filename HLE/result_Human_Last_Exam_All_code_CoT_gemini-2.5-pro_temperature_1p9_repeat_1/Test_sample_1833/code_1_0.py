import math

def calculate_tension():
    """
    This function calculates the tension in the string for the given physics problem.
    """
    # Step 1: Define the given parameters
    m = 1.0  # mass of the ring in kg
    M = 1.0  # mass of the object in kg
    g = 9.8  # acceleration due to gravity in m/s^2
    theta_deg = 60.0  # angle in degrees

    # Convert angle to radians for use in trigonometric functions
    theta_rad = math.radians(theta_deg)
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)
    cos_sq_theta = cos_theta**2

    # Step 2: Calculate the term L*omega^2 using conservation of energy.
    # L*omega^2 represents the centripetal acceleration of M relative to the ring m, multiplied by L.
    # The formula is derived from energy and momentum conservation:
    # L*omega^2 = (2 * g * sin(theta) * (m + M)) / (m + M * cos(theta)^2)
    numerator_L_omega_sq = 2 * g * sin_theta * (m + M)
    denominator_L_omega_sq = m + M * cos_sq_theta
    L_omega_sq = numerator_L_omega_sq / denominator_L_omega_sq

    # Step 3: Calculate the acceleration of the ring, a_m.
    # This is found by analyzing the forces and time-differentiating the momentum conservation equation.
    # The formula is: a_m = (M * cos(theta) * (g * sin(theta) + L*omega^2)) / (m + M * cos(theta)^2)
    numerator_a_m = M * cos_theta * (g * sin_theta + L_omega_sq)
    denominator_a_m = m + M * cos_sq_theta
    a_m = numerator_a_m / denominator_a_m

    # Step 4: Calculate the tension T using the radial force equation on mass M.
    # The equation for tension is: T = T_gravity + T_centripetal + T_ring_acc
    # which corresponds to the formula: T = M*g*sin(theta) + M*L*omega^2 + M*a_m*cos(theta)

    term1_gravity = M * g * sin_theta
    term2_centripetal = M * L_omega_sq
    term3_ring_acc = M * a_m * cos_theta

    tension = term1_gravity + term2_centripetal + term3_ring_acc

    # Step 5: Output the results as requested.
    print(f"The equation for the tension is T = M*g*sin(theta) + M*L*omega^2 + M*a_m*cos(theta)")
    print(f"For theta = {theta_deg} degrees, the components of tension are:")
    print(f"T = {round(term1_gravity, 2)} N (from gravity) + {round(term2_centripetal, 2)} N (from relative motion) + {round(term3_ring_acc, 2)} N (from ring's acceleration)")
    print(f"\nThe calculated total tension in the string is {round(tension, 2)} N.")

if __name__ == '__main__':
    calculate_tension()