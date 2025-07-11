import math

def solve_drawbridge_velocity():
    """
    Calculates the vertical velocity of a drawbridge edge based on given parameters.
    """
    # Given parameters from the problem
    L = 40  # Length of the drawbridge in meters
    y = 10  # Vertical height of the bridge edge in meters
    
    # The relationship between vertical height y, length L, and angle θ is:
    # y = L * cos(θ)
    
    # The vertical velocity is dy/dt. By differentiating, we get:
    # dy/dt = -L * sin(θ) * dθ/dt
    
    print("Step 1: Find the angle θ when the bridge is 10 m high.")
    # From y = L * cos(θ), we can find cos(θ)
    cos_theta = y / L
    # From sin²(θ) + cos²(θ) = 1, we find sin(θ).
    # Since the bridge is being raised, 0 < θ < π/2, so sin(θ) is positive.
    sin_theta = math.sqrt(1 - cos_theta**2)
    
    print(f"The length of the bridge L = {L} m.")
    print(f"The height of the edge y = {y} m.")
    print(f"From y = L * cos(θ), we find cos(θ) = {y}/{L} = {cos_theta}.")
    print(f"From sin²(θ) + cos²(θ) = 1, we find sin(θ) = sqrt(1 - {cos_theta}²) = {sin_theta:.4f}.")
    print("-" * 20)
    
    print("Step 2: Use the given rate of change for θ, dθ/dt.")
    # The problem provides dθ/dt = -(3π/10) / cos(π/12)
    # Let's calculate this value.
    cos_pi_12 = math.cos(math.pi / 12)
    dtheta_dt_val = -(3 * math.pi / 10) / cos_pi_12
    
    print("The given rate dθ/dt = -(3 * π / 10) / cos(π / 12).")
    print(f"The value is dθ/dt = -(3 * {math.pi:.4f} / 10) / {cos_pi_12:.4f} = {dtheta_dt_val:.4f} rad/s.")
    print("-" * 20)

    print("Step 3: Calculate the vertical velocity dy/dt.")
    # The formula for vertical velocity is dy/dt = -L * sin(θ) * dθ/dt
    # Now we substitute the values we found.
    dy_dt = -L * sin_theta * dtheta_dt_val
    
    print("The equation for the vertical velocity is: dy/dt = -L * sin(θ) * dθ/dt")
    print("Substituting the numbers into the final equation:")
    print(f"dy/dt = -({L}) * ({sin_theta:.4f}) * ({dtheta_dt_val:.4f})")
    print(f"The final calculated vertical velocity is: {dy_dt:.4f} m/s.")

solve_drawbridge_velocity()