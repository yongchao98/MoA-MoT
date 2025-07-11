import math

def solve_drawbridge_velocity():
    """
    Calculates the vertical velocity of the edge of a drawbridge.
    """
    # --- Given values ---
    L = 40.0  # Length of the drawbridge in meters
    y = 10.0  # Height of the bridge end in meters

    # --- Step 1: Find geometric properties at the given moment ---
    # The relationship is y = L * cos(theta).
    # We find cos(theta) when y = 10.
    cos_theta = y / L
    
    # Find sin(theta) using the identity sin^2(theta) + cos^2(theta) = 1.
    # We take the positive root as 0 <= theta <= pi/2.
    sin_theta = math.sqrt(1 - cos_theta**2)

    # --- Step 2: Use the given rate of change for the angle theta ---
    # d(theta)/dt = -(3*pi/10) / cos(pi/12)
    d_theta_dt = -(3 * math.pi / 10) / math.cos(math.pi / 12)

    # --- Step 3: Apply the related rates formula ---
    # The velocity formula is dy/dt = -L * sin(theta) * d(theta)/dt
    dy_dt = -L * sin_theta * d_theta_dt

    # --- Step 4: Output the equation with values and the final answer ---
    print("The vertical velocity (dy/dt) is calculated using the formula:")
    print("dy/dt = -L * sin(\u03B8) * d\u03B8/dt\n")
    print("At the moment when the height y = 10 m:")
    print(f"  L = {L}")
    print(f"  sin(\u03B8) = \u221A(1 - (1/4)\u00b2) = {sin_theta:.4f}")
    print(f"  d\u03B8/dt = -(3\u03C0/10) / cos(\u03C0/12) = {d_theta_dt:.4f}\n")
    
    print("Plugging the numbers into the equation:")
    print(f"dy/dt = -({L}) * ({sin_theta:.4f}) * ({d_theta_dt:.4f})")
    
    print("\nThe final calculated vertical velocity is:")
    print(f"{dy_dt:.4f} m/s")

solve_drawbridge_velocity()