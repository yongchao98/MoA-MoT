import math

def solve_tension():
    """
    Calculates the tension in the string for the given physics problem.
    """
    # Given parameters
    m = 1.0  # mass of the ring in kg
    M = 1.0  # mass of the object in kg
    g = 9.8  # acceleration due to gravity in m/s^2
    theta_deg = 60.0  # angle in degrees

    # Convert angle to radians for trigonometric functions
    theta_rad = math.radians(theta_deg)

    # Pre-calculate trigonometric values
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)
    cot_theta = 1.0 / math.tan(theta_rad)
    sin2_theta = sin_theta**2
    cos2_theta = cos_theta**2

    # Since m = M, the equations simplify.
    # Step 1: Calculate the term L*\dot{θ}^2 from energy conservation.
    # The term represents the contribution from centripetal-like motion.
    # L * (dθ/dt)^2 = (4 * g * sin(θ)) / (2 - sin(θ)^2)
    L_thetadot2 = (4 * g * sin_theta) / (2 - sin2_theta)
    
    # Step 2: Calculate the term L*\ddot{θ} from the equation of motion.
    # This term represents the contribution from angular acceleration.
    # L * (d^2θ/dt^2) = (g*cos(θ)*(2*sin(θ)^2 + 4)) / ((1 + cos(θ)^2)*(2-sin(θ)^2))
    L_thetaddot_num = g * cos_theta * (2 * sin2_theta + 4)
    L_thetaddot_den = (1 + cos2_theta) * (2 - sin2_theta)
    L_thetaddot = L_thetaddot_num / L_thetaddot_den
    
    # Step 3: Calculate the three components of the tension formula:
    # T = M*g/sin(θ) - M*(L*\ddot{θ})*cot(θ) + M*(L*\dot{θ}^2)
    
    # Component from gravity projected onto the vertical direction
    term_gravity = M * g / sin_theta
    
    # Component from the tangential acceleration part
    term_accel = M * L_thetaddot * cot_theta
    
    # Component from the radial acceleration part
    term_velocity = M * L_thetadot2

    # Calculate the final tension
    T = term_gravity - term_accel + term_velocity

    # Print the equation with the calculated values of each term
    print(f"The equation for tension is: T = Mg/sin(θ) - M(L\u00d7\u00a8)\u00d7cot(θ) + M(L\u00d7\u02d9θ\u00b2)")
    print(f"Calculating each term with M=1 kg, g=9.8 m/s^2, θ=60°:")
    print(f"First term (gravity component): {term_gravity:.2f} N")
    print(f"Second term (angular acceleration component): {term_accel:.2f} N")
    print(f"Third term (angular velocity component): {term_velocity:.2f} N")
    print(f"Tension T = {term_gravity:.2f} - {term_accel:.2f} + {term_velocity:.2f}")

    # Print the final result
    print(f"\nThe tension in the string is: {T:.2f} Newtons.")

    return T

# Run the calculation and store the result
tension_value = solve_tension()

# The final answer in the required format
# <<<28.52>>>
print(f"\n<<<{tension_value:.2f}>>>")