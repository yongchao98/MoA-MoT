import math

def solve_tension():
    """
    Calculates the tension in the string for the given physics problem.
    """
    # 1. Define given constants and initial conditions
    m = 1.0  # mass of the ring in kg
    M = 1.0  # mass of the object in kg
    g = 9.8  # acceleration due to gravity in m/s^2
    theta_deg = 60.0
    theta_rad = math.radians(theta_deg)

    # Pre-calculate trigonometric values for efficiency
    cos_theta = math.cos(theta_rad)
    sin_theta = math.sin(theta_rad)

    # 2. Calculate L*ω² (L*omega_squared) using energy and momentum conservation.
    # The length L cancels out from all equations.
    # Lω² = 2*g*sin(θ)*(m+M) / (m + M*cos²θ)
    l_omega_sq = (2 * g * sin_theta * (m + M)) / (m + M * cos_theta**2)

    # 3. Calculate L*α (L*alpha) using the tangential force equation.
    # Lα = [(m+M)/(m+Mcos²θ)] * [g*cosθ + (M*Lω²/(m+M))*cosθ*sinθ]
    term1_alpha = (m + M) / (m + M * cos_theta**2)
    term2_alpha = g * cos_theta + (M * l_omega_sq / (m + M)) * cos_theta * sin_theta
    l_alpha = term1_alpha * term2_alpha

    # 4. Calculate the acceleration of the ring, a_m.
    # a_m = (M/(m+M)) * [Lω²*cosθ + Lα*sinθ]
    am = (M / (m + M)) * (l_omega_sq * cos_theta + l_alpha * sin_theta)

    # 5. Calculate the tension T using the radial force equation.
    # T = M*(g*sinθ - a_m*cosθ + Lω²)
    T = M * (g * sin_theta - am * cos_theta + l_omega_sq)

    # 6. Print the results, showing the final calculation as requested.
    print("This script calculates the tension in the string based on the derived physical formulas.\n")
    print(f"Given values: m = {m} kg, M = {M} kg, g = {g} m/s^2, θ = {theta_deg}°\n")
    print("Intermediate calculated values:")
    print(f"  Lω² (related to kinetic energy): {l_omega_sq:.4f} m/s²")
    print(f"  a_m (acceleration of the ring): {am:.4f} m/s²\n")
    
    print("Final tension calculation:")
    print(f"T = M * (g * sin(θ) - a_m * cos(θ) + Lω²)")
    print(f"T = {M:.1f} * ({g:.1f} * {sin_theta:.4f} - {am:.4f} * {cos_theta:.4f} + {l_omega_sq:.4f})")
    
    g_sin_theta = g * sin_theta
    am_cos_theta = am * cos_theta
    print(f"T = {M:.1f} * ({g_sin_theta:.4f} - {am_cos_theta:.4f} + {l_omega_sq:.4f})")
    
    result_in_parentheses = g_sin_theta - am_cos_theta + l_omega_sq
    print(f"T = {M:.1f} * ({result_in_parentheses:.4f})")
    
    print(f"\nThe tension in the string is {T:.2f} Newtons.")

solve_tension()