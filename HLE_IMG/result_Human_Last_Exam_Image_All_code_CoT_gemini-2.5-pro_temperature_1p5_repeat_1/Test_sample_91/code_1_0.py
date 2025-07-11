import math

def solve_drawbridge_velocity():
    """
    Calculates the vertical velocity of a drawbridge edge at a specific height.
    """
    # --- Given values from the problem statement ---
    L = 40  # Length of the drawbridge in meters
    y_target = 10  # Target vertical height in meters

    # --- Step 1: Find cos(theta) and sin(theta) at the target height ---
    # The vertical height 'y' is given by y = L * cos(theta).
    cos_theta = y_target / L
    
    # We find sin(theta) using the identity sin^2(theta) + cos^2(theta) = 1.
    # As the bridge is raised, theta is in (0, pi/2), so sin(theta) is positive.
    sin_theta = math.sqrt(1 - cos_theta**2)

    # --- Step 2: Use the given rate of change of theta ---
    # The rate of change of theta, d(theta)/dt, is given.
    # d(theta)/dt = - (3*pi/10) / cos(pi/12)
    dtheta_dt_val = - (3 * math.pi / 10) / math.cos(math.pi / 12)

    # --- Step 3: Find the vertical velocity dy/dt ---
    # The vertical velocity, dy/dt, is the time derivative of y = L * cos(theta).
    # Using the chain rule: dy/dt = -L * sin(theta) * d(theta)/dt
    vertical_velocity = -L * sin_theta * dtheta_dt_val

    # --- Step 4: Print the calculation step-by-step ---
    print("The final equation for the vertical velocity (dy/dt) is derived from y = L * cos(theta).")
    print("dy/dt = -L * sin(theta) * d(theta)/dt")
    print("\nSubstituting the specific values:")
    
    # Print each number in the final equation as requested.
    print(f"L = {L}")
    print(f"sin(theta) for y={y_target} is {sin_theta:.4f}")
    print(f"d(theta)/dt is given as a constant value {dtheta_dt_val:.4f} rad/min")
    print("\nFinal calculation:")
    print(f"dy/dt = -({L}) * ({sin_theta:.4f}) * ({dtheta_dt_val:.4f})")
    print(f"dy/dt = {vertical_velocity:.2f} m/min")

solve_drawbridge_velocity()
<<<37.78>>>