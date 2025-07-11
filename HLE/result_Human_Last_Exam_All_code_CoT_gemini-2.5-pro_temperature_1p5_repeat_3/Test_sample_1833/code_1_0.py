import math

def calculate_tension():
    """
    Calculates the tension in the string for the given physics problem.
    """
    # --- 1. Define constants and initial values ---
    m = 1.0  # kg, mass of the ring
    M = 1.0  # kg, mass of the object
    g = 9.8  # m/s^2, acceleration due to gravity
    theta_deg = 60.0 # degrees

    # --- 2. Convert angle to radians and calculate sin/cos ---
    theta_rad = math.radians(theta_deg)
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)

    # --- 3. Calculate L*(d(theta)/dt)^2 using energy conservation ---
    # This term is related to the kinetic energy of the system.
    # Formula derived from Conservation of Energy:
    # L*(d(theta)/dt)^2 = 2 * g * sin(theta) / [m * sin(theta)^2 / (m + M) + cos(theta)^2]
    l_theta_dot_sq_num = 2 * g * sin_theta
    l_theta_dot_sq_den = (m * sin_theta**2) / (m + M) + cos_theta**2
    l_theta_dot_sq = l_theta_dot_sq_num / l_theta_dot_sq_den

    # --- 4. Calculate Tension (T) using force analysis ---
    # Formula derived from Newton's Second Law in a non-inertial frame:
    # T = (m*M / (m + M*cos(theta)^2)) * (g*sin(theta) + L*(d(theta)/dt)^2)
    prefactor_num = m * M
    prefactor_den = m + M * cos_theta**2
    prefactor = prefactor_num / prefactor_den
    paren_term = g * sin_theta + l_theta_dot_sq
    T = prefactor * paren_term

    # --- 5. Print the step-by-step calculation and the final answer ---
    print("Step-by-step calculation for the tension in the string:")
    print("-" * 60)
    print(f"Given values: m = {m:.1f} kg, M = {M:.1f} kg, g = {g:.1f} m/s^2, theta = {theta_deg:.1f} degrees")
    print(f"(sin(60) = {sin_theta:.3f}, cos(60) = {cos_theta:.3f})")
    print("")

    # Step 1: Calculation of the speed term
    print("1. First, find the term L*(d(theta)/dt)^2 from energy conservation:")
    print(f"   L*(d(theta)/dt)^2 = (2 * {g:.1f} * {sin_theta:.3f}) / [({m:.1f} * {sin_theta:.3f}^2)/({m:.1f}+{M:.1f}) + {cos_theta:.3f}^2]")
    print(f"                     = {l_theta_dot_sq_num:.4f} / {l_theta_dot_sq_den:.4f}")
    print(f"                     = {l_theta_dot_sq:.2f}")
    print("")

    # Step 2: Calculation of Tension
    print("2. Second, calculate the tension T using the force equation:")
    print(f"   T = ({m:.1f}*{M:.1f} / ({m:.1f} + {M:.1f}*{cos_theta:.3f}^2)) * ({g:.1f}*{sin_theta:.3f} + {l_theta_dot_sq:.2f})")
    print(f"     = ({prefactor:.3f}) * ({paren_term:.3f})")
    print(f"\nThe final tension in the string is: {T:.2f} N")
    print("-" * 60)

if __name__ == '__main__':
    calculate_tension()