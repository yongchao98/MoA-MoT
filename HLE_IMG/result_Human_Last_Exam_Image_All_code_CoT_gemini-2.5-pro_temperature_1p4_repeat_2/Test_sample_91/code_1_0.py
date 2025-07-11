import math

def solve_drawbridge_velocity():
    """
    Calculates the vertical velocity of a drawbridge edge at a specific height.
    """
    # --- Given values ---
    L = 40.0  # Length of the drawbridge in meters
    y = 10.0   # Vertical height of the edge in meters
    
    # --- Plan ---
    # 1. Find sin(theta) from y = L * cos(theta)
    # 2. Calculate the given constant d(theta)/dt
    # 3. Calculate dy/dt = -L * sin(theta) * d(theta)/dt

    # --- Step 1: Find sin(theta) ---
    # From y = L * cos(theta), we have cos(theta) = y / L
    cos_theta = y / L
    
    # Using sin^2(theta) + cos^2(theta) = 1, we get sin(theta) = sqrt(1 - cos^2(theta))
    # Since the bridge is being raised from horizontal, 0 < theta < pi/2, so sin(theta) is positive.
    sin_theta = math.sqrt(1 - cos_theta**2)

    # --- Step 2: Calculate d(theta)/dt ---
    # The rate of change of theta is given as a constant value.
    # d(theta)/dt = -(3 * pi / 10) / cos(pi / 12)
    dtheta_dt = -(3 * math.pi / 10) / math.cos(math.pi / 12)

    # --- Step 3: Calculate dy/dt ---
    # The vertical velocity is the derivative of y with respect to t:
    # dy/dt = d/dt (L * cos(theta)) = -L * sin(theta) * d(theta)/dt
    dy_dt = -L * sin_theta * dtheta_dt

    # --- Output the results as requested ---
    print("The formula for vertical velocity is: dy/dt = -L * sin(theta) * d(theta)/dt")
    print("\nWe are given the following values:")
    print(f"Length of the bridge, L = {L} m")
    
    print("\nAt the specific moment when the vertical height is y = 10 m, we find:")
    print(f"sin(theta) = {sin_theta:.4f}")
    
    print("\nThe given rate of change for the angle is a constant:")
    print(f"d(theta)/dt = {dtheta_dt:.4f} rad/s")

    print("\nSubstituting these values into the velocity equation:")
    print(f"dy/dt = -({L}) * ({sin_theta:.4f}) * ({dtheta_dt:.4f})")
    
    print("\nThe calculated vertical velocity is:")
    print(f"dy/dt = {dy_dt:.4f} m/s")

solve_drawbridge_velocity()