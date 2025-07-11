import math

def solve_drawbridge_velocity():
    """
    Calculates the vertical velocity of a drawbridge at a specific height.
    """
    # Given parameters from the problem
    L = 40  # Length of the drawbridge in meters
    h = 10  # Height of the bridge's edge in meters
    
    # The problem provides the rate of change of the angle theta
    # dθ/dt = - (3 * π / 10) / cos(π / 12)
    # We will calculate this value.
    d_theta_dt = -(3 * math.pi / 10) / math.cos(math.pi / 12)

    # --- Step 1: Find cos(θ) and sin(θ) at the given instant ---
    # The geometric relationship is h = L * cos(θ)
    # So, cos(θ) = h / L
    cos_theta = h / L
    
    # Using the identity sin²(θ) + cos²(θ) = 1, we find sin(θ).
    # sin(θ) = sqrt(1 - cos²(θ)). Since the bridge is rising from horizontal,
    # θ is between 0 and π/2, so sin(θ) is positive.
    sin_theta = math.sqrt(1 - cos_theta**2)

    # --- Step 2: Calculate the vertical velocity dh/dt ---
    # Differentiating h = L * cos(θ) with respect to time t gives:
    # dh/dt = -L * sin(θ) * dθ/dt
    dh_dt = -L * sin_theta * d_theta_dt

    # --- Step 3: Output the results as requested ---
    print("The equation for the vertical velocity is: dh/dt = -L * sin(θ) * dθ/dt\n")
    
    print("We are given or can calculate the following values:")
    print(f"L = {L}")
    print(f"h = {h}")
    # Value of sin(θ) when h=10. sin(θ) = sqrt(1 - (10/40)^2) = sqrt(1 - 1/16) = sqrt(15)/4
    print(f"sin(θ) = sqrt(15)/4 ≈ {sin_theta:.4f}")
    print(f"dθ/dt = -(3π/10)/cos(π/12) ≈ {d_theta_dt:.4f}\n")

    print("Plugging these values into the velocity equation:")
    print(f"dh/dt = -({L}) * ({sin_theta:.4f}) * ({d_theta_dt:.4f})")

    print(f"\nThe final vertical velocity of the moving edge of the bridge is: {dh_dt:.2f} m/s")

solve_drawbridge_velocity()