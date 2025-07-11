import math

def solve_drawbridge_velocity():
    """
    Calculates the vertical velocity of the drawbridge edge.
    """
    # Given parameters from the problem
    L = 40  # Length of the drawbridge in meters
    y = 10  # Vertical height of the edge of the bridge in meters

    # Given rate of change of theta (dθ/dt).
    # Based on the problem context ("revolutions per minute"), we assume the time unit is minutes.
    d_theta_dt = -(3 * math.pi) / (10 * math.cos(math.pi / 12))

    # The relationship between vertical height y and angle θ is y = L * cos(θ).
    # The vertical velocity dy/dt is its derivative: dy/dt = -L * sin(θ) * dθ/dt.

    # Step 1: Find cos(θ) at the moment y = 10 m.
    cos_theta = y / L

    # Step 2: Find sin(θ) from cos(θ).
    # sin(θ) is positive as the bridge is being raised from horizontal (0 < θ < π/2).
    sin_theta = math.sqrt(1 - cos_theta**2)

    # Step 3: Calculate the vertical velocity dy/dt.
    dy_dt = -L * sin_theta * d_theta_dt

    # Print the steps and the final equation with values
    print("To find the vertical velocity (dy/dt), we use the formula: dy/dt = -L * sin(θ) * dθ/dt")
    print("\nGiven values and calculated intermediates:")
    print(f"L = {L}")
    print(f"When y = {y}, cos(θ) = {y}/{L} = {cos_theta}")
    print(f"sin(θ) = sqrt(1 - {cos_theta}^2) = {sin_theta:.4f}")
    print(f"The given dθ/dt = {d_theta_dt:.4f} rad/min")
    
    print("\nSubstituting these values into the velocity equation:")
    print(f"dy/dt = -({L}) * ({sin_theta:.4f}) * ({d_theta_dt:.4f})")
    
    # Final result
    print(f"\nThe vertical velocity of the moving edge of the bridge is {dy_dt:.4f} m/min.")
    print(f"<<<{dy_dt:.4f}>>>")

solve_drawbridge_velocity()