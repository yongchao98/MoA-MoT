import math

def solve_drawbridge_velocity():
    """
    Calculates the vertical velocity of a drawbridge at a specific moment.
    """
    # --- Step 1: Define the given parameters ---
    L = 40.0  # Length of the drawbridge in meters
    h = 10.0  # Height of the bridge edge at the specific moment in meters

    # The rate of change of θ is given. Note: pi/12 radians is 15 degrees.
    dtheta_dt = -(3 * math.pi / 10) / math.cos(math.pi / 12)

    # --- Step 2: Find sin(θ) at the given instant (h = 10) ---
    # The geometric relationship is h = L * cos(θ), so cos(θ) = h / L
    cos_theta = h / L

    # Using the identity sin²(θ) + cos²(θ) = 1, we find sin(θ).
    # Since the bridge is rising from horizontal, θ is between 0 and π/2, so sin(θ) is positive.
    sin_theta = math.sqrt(1 - cos_theta**2)

    # --- Step 3: Calculate the vertical velocity (dh/dt) ---
    # The relationship between the rates is found by differentiating h = L * cos(θ)
    # This gives: dh/dt = -L * sin(θ) * dθ/dt
    dh_dt = -L * sin_theta * dtheta_dt

    # --- Step 4: Print the results and the final equation ---
    print("To find the vertical velocity (dh/dt), we use the formula: dh/dt = -L * sin(θ) * dθ/dt")
    print("\nAt the moment when h = 10 m, we have:")
    print(f"L = {L:.1f}")
    print(f"sin(θ) = {sin_theta:.6f}")
    print(f"dθ/dt = {dtheta_dt:.6f}")
    
    print("\nPlugging these numbers into the final equation:")
    # The problem requires printing each number in the final equation.
    print(f"dh/dt = -({L:.1f}) * ({sin_theta:.6f}) * ({dtheta_dt:.6f})")

    print(f"\nThe vertical velocity of the moving edge is {dh_dt:.4f} m/s.")

solve_drawbridge_velocity()