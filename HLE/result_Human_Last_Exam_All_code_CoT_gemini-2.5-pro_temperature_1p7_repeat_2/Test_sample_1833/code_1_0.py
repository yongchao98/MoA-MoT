import math

def calculate_tension():
    """
    Calculates the tension in the string for the given physics problem.
    """
    # --- Given parameters ---
    m = 1.0  # kg, mass of the ring
    M = 1.0  # kg, mass of the object
    g = 9.8  # m/s^2, acceleration due to gravity
    theta_deg = 60.0  # degrees

    # --- Convert angle to radians for math functions ---
    theta_rad = math.radians(theta_deg)

    # --- Calculate intermediate trigonometric values for clarity ---
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)
    cos_sq_theta = cos_theta**2

    # --- Print explanation and formulas ---
    print("This problem is solved by combining conservation of energy and Newton's second law.")
    
    print("\nStep 1: Use conservation of energy to find the squared angular velocity term (L*omega^2).")
    print("The derived relationship is:")
    print("L*omega^2 = (2 * g * sin(theta) * (m + M)) / (m + M * cos(theta)^2)")
    print(f"L*omega^2 = (2 * {g:.1f} * sin({theta_deg}) * ({m:.1f} + {M:.1f})) / ({m:.1f} + {M:.1f} * cos({theta_deg})^2)")
    
    # --- Step 1 Calculation ---
    L_omega_sq_num = 2 * g * sin_theta * (m + M)
    L_omega_sq_den = m + M * cos_sq_theta
    L_omega_sq = L_omega_sq_num / L_omega_sq_den
    
    print(f"Result of Step 1: L*omega^2 = {L_omega_sq:.4f}")

    print("\nStep 2: Use the equations of motion to derive and calculate the tension T.")
    print("The derived formula for tension is:")
    print("T = (M * (g * sin(theta) + L*omega^2)) / (1 + (M/m) * cos(theta)^2)")
    print(f"T = ({M:.1f} * ({g:.1f} * sin({theta_deg}) + {L_omega_sq:.4f})) / (1 + ({M:.1f}/{m:.1f}) * cos({theta_deg})^2)")

    # --- Step 2 Calculation ---
    T_num = M * (g * sin_theta + L_omega_sq)
    T_den = 1 + (M / m) * cos_sq_theta
    tension = T_num / T_den
    
    # Round the final result to two decimal places
    tension_rounded = round(tension, 2)
    
    print(f"\nFinal calculation gives T = {tension:.4f} Newtons.")
    print(f"The tension in the string, rounded to two decimal places, is: {tension_rounded} N.")


if __name__ == "__main__":
    calculate_tension()
    print("\n<<<28.52>>>")
