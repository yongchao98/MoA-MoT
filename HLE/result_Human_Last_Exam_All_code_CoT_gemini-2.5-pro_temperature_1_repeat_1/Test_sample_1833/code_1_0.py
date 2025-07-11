import math

def calculate_tension():
    """
    Calculates the tension in the string for the given physics problem.

    The solution involves applying conservation of energy and conservation of 
    horizontal momentum to find the velocities and accelerations of the system 
    at the specified angle. Then, Newton's second law is used in the 
    non-inertial frame of the ring to find the tension.
    """
    # Given values
    m = 1.0  # mass of the ring in kg
    M = 1.0  # mass of the object in kg
    g = 9.8  # acceleration due to gravity in m/s^2
    theta_deg = 60.0 # angle in degrees

    # Convert angle to radians for trigonometric functions
    theta_rad = math.radians(theta_deg)

    # Pre-calculate sine and cosine of the angle
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)

    # 1. Use conservation of energy and momentum to find (L * (d(theta)/dt)^2).
    # This term represents the square of the object's speed relative to the ring, divided by L.
    # L * (d(theta)/dt)^2 = (2 * g * sin(theta) * (m + M)) / (m + M * cos(theta)^2)
    # This is derived from combining the energy and horizontal momentum conservation equations.
    # The length L cancels out from the final expressions.
    denominator = m + M * cos_theta**2
    L_theta_dot_sq = (2 * g * sin_theta * (m + M)) / denominator

    # 2. Find the horizontal acceleration of the ring (a_m).
    # This is found by analyzing the forces and accelerations in the system.
    # a_m = (M * cos(theta) * (L_theta_dot_sq - g * sin(theta))) / (m + M * cos(theta)^2)
    a_m = (M * cos_theta * (L_theta_dot_sq - g * sin_theta)) / denominator

    # 3. Calculate the tension (T) by analyzing forces on mass M in the ring's non-inertial frame.
    # The tension T balances three effective forces in the radial direction:
    # - The component of gravity along the string: M*g*sin(theta)
    # - The component of the inertial force (-M*a_m) along the string: M*a_m*cos(theta)
    # - The centripetal force required for circular motion relative to the ring: M*L*(d(theta)/dt)^2
    
    # Calculate each component of the tension equation T = term1 + term2 + term3
    gravitational_component = M * g * sin_theta
    inertial_component = M * a_m * cos_theta
    centripetal_component = M * L_theta_dot_sq
    
    # Total tension is the sum of these components
    T = gravitational_component + inertial_component + centripetal_component

    # Print the final equation with numerical values
    print("The equation for tension T is the sum of the radial components of gravity, inertial force, and the centripetal force:")
    print(f"T = {gravitational_component:.2f} N (gravity) + {inertial_component:.2f} N (inertial) + {centripetal_component:.2f} N (centripetal)")
    
    # Print the final result
    print(f"The final equation is: {gravitational_component:.2f} + {inertial_component:.2f} + {centripetal_component:.2f} = {T:.2f}")
    print(f"The tension in the string is {T:.2f} Newtons.")


if __name__ == "__main__":
    calculate_tension()