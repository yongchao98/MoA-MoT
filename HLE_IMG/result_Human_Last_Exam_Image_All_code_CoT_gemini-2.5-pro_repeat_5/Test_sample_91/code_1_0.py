import math

def solve_drawbridge_velocity():
    """
    Calculates the vertical velocity of the edge of a drawbridge.
    """
    # --- Given values ---
    L = 40  # Length of the drawbridge in meters
    y = 10  # Vertical height of the edge of the bridge in meters

    # --- Step 1: Establish the geometric relationship and find sin(θ) ---
    # The relationship between height y, length L, and angle θ is:
    # y = L * cos(θ)
    # We can find cos(θ) when y = 10.
    cos_theta = y / L
    # Using the identity sin²(θ) + cos²(θ) = 1, we find sin(θ).
    # Since the bridge is being raised from horizontal, θ is between 0 and π/2, so sin(θ) is positive.
    sin_theta = math.sqrt(1 - cos_theta**2)

    # --- Step 2: Use the given rate of change of θ ---
    # The rate of change of theta is given as dθ/dt = -(3π/10) / cos(π/12).
    # Let's calculate this value.
    cos_pi_12 = math.cos(math.pi / 12)
    dtheta_dt_val = -(3 * math.pi / 10) / cos_pi_12

    # --- Step 3: Find the vertical velocity dy/dt ---
    # Differentiate the position equation y = L * cos(θ) with respect to time t.
    # The result is dy/dt = -L * sin(θ) * dθ/dt (using the chain rule).
    # Now, substitute the known values to find the velocity.
    dy_dt = -L * sin_theta * dtheta_dt_val
    
    # --- Step 4: Print the explanation and the final equation with numbers ---
    print("The relationship between the vertical height (y), length (L), and angle (θ) is:")
    print(f"y = L * cos(θ)")
    print("\nDifferentiating with respect to time t gives the velocity relationship:")
    print("dy/dt = -L * sin(θ) * dθ/dt")

    print("\nAt the instant when the height y = 10 m:")
    print(f"The length of the bridge L = {L} m.")
    print(f"cos(θ) = y / L = {y}/{L} = {cos_theta}")
    print(f"sin(θ) = sqrt(1 - ({cos_theta})²) = {sin_theta:.5f}")
    
    print("\nThe given rate of change of the angle is:")
    # We are given dθ/dt = -(3π/10) / cos(π/12)
    print(f"dθ/dt = -(3 * π / 10) / cos(π / 12) ≈ {dtheta_dt_val:.5f} rad/s")

    print("\nSubstituting these values into the velocity equation:")
    print(f"dy/dt = -({L}) * ({sin_theta:.5f}) * ({dtheta_dt_val:.5f})")
    
    print("\nTherefore, the calculated vertical velocity of the bridge's edge is:")
    print(f"dy/dt = {dy_dt:.4f} m/s")

    # Final answer for the system
    return dy_dt

# Run the solver
final_velocity = solve_drawbridge_velocity()
# <<<37.7884>>>