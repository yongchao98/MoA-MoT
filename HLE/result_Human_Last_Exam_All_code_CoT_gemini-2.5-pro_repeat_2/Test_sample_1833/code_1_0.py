import math

def calculate_tension():
    """
    Calculates the tension in the string for the given physics problem.
    """
    # Given parameters
    m = 1.0  # mass of the ring in kg
    M = 1.0  # mass of the object in kg
    g = 9.8  # acceleration due to gravity in m/s^2
    theta_deg = 60.0  # angle in degrees

    # Convert angle to radians
    theta_rad = math.radians(theta_deg)

    # Trigonometric values
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)
    cot_theta = 1.0 / math.tan(theta_rad)

    # Step 1: Calculate L*(d(theta)/dt)^2 from energy conservation
    # (1/2)*K*L^2*(d(theta)/dt)^2 = M*g*L*sin(theta)
    # where K = (M/(m+M))*(m + M*cos(theta)^2)
    # This simplifies to: L*(d(theta)/dt)^2 = (2*g*(m+M)*sin(theta)) / (m + M*cos(theta)^2)
    L_theta_dot_sq = (2 * g * (m + M) * sin_theta) / (m + M * cos_theta**2)

    # Step 2: Calculate L*(d^2(theta)/dt^2) from the equation of motion
    # L*(d^2(theta)/dt^2) = (M*L*(d(theta)/dt)^2*sin(theta)*cos(theta) + g*(m+M)*cos(theta)) / (m + M*cos(theta)^2)
    L_theta_ddot_numerator = M * L_theta_dot_sq * sin_theta * cos_theta + g * (m + M) * cos_theta
    L_theta_ddot_denominator = m + M * cos_theta**2
    L_theta_ddot = L_theta_ddot_numerator / L_theta_ddot_denominator
    
    # Step 3: Calculate Tension T using the force equation in the y-direction for mass M
    # T*sin(theta) - M*g = M*(-L*(d^2(theta)/dt^2)*cos(theta) + L*(d(theta)/dt)^2*sin(theta))
    # Rearranging for T:
    # T = M*g/sin(theta) + M*L*(d(theta)/dt)^2 - M*L*(d^2(theta)/dt^2)*cot(theta)
    
    term1 = M * g / sin_theta
    term2 = M * L_theta_dot_sq
    term3 = M * L_theta_ddot * cot_theta
    
    tension = term1 + term2 - term3

    print("The final equation for tension is T = A + B - C, where:")
    print(f"A (gravity component) = {M:.1f} * {g:.1f} / sin({theta_deg:.0f}) = {term1:.2f} N")
    print(f"B (centripetal-related component) = {M:.1f} * (L*(dθ/dt)²) = {term2:.2f} N")
    print(f"C (tangential acceleration component) = {M:.1f} * (L*(d²θ/dt²)) * cot({theta_deg:.0f}) = {term3:.2f} N")
    print("\nFinal Calculation:")
    print(f"T = {term1:.2f} + {term2:.2f} - {term3:.2f}")
    print(f"The tension in the string is: {tension:.2f} N")

calculate_tension()