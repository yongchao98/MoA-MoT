import math

def solve_drawbridge_velocity():
    """
    Calculates the vertical velocity of a drawbridge edge at a specific height.
    """
    # Given parameters from the problem
    L = 40.0  # Length of the drawbridge in meters
    y = 10.0  # Vertical height of the bridge edge in meters
    
    # The rate of change of theta is given as a constant value.
    # The time unit is assumed to be in minutes, to be consistent with the 
    # extraneous information about wheel rotation.
    dtheta_dt_val = -(3 * math.pi / 10) / math.cos(math.pi / 12)

    # Step 1: Find sin(theta) at the moment when y = 10 m.
    # The relationship is y = L * cos(theta), so cos(theta) = y / L.
    cos_theta = y / L
    # Using the identity sin^2(theta) + cos^2(theta) = 1.
    # sin(theta) must be positive as the bridge is rising.
    sin_theta = math.sqrt(1 - cos_theta**2)

    # Step 2: The vertical velocity is dy/dt, found by differentiating y = L * cos(theta).
    # dy/dt = -L * sin(theta) * d(theta)/dt
    vertical_velocity = -L * sin_theta * dtheta_dt_val
    
    # Step 3: Print the explanation and the final equation with all numbers.
    print("The equation for the vertical velocity is: dy/dt = -L * sin(theta) * d(theta)/dt")
    print("\nWe are given the following values:")
    print(f"Length of the bridge, L = {L} m")
    print(f"Rate of change of the angle, d(theta)/dt = -(3*pi/10)/cos(pi/12) ≈ {dtheta_dt_val:.4f} rad/min")

    print("\nAt the moment when the vertical height y = 10 m:")
    print(f"cos(theta) = y / L = {y} / {L} = {cos_theta}")
    print(f"sin(theta) = sqrt(1 - cos^2(theta)) = sqrt(1 - {cos_theta**2}) ≈ {sin_theta:.4f}")

    print("\nSubstituting these numbers into the velocity equation:")
    # The prompt requires printing each number in the final equation.
    print(f"dy/dt = -({L}) * ({sin_theta:.4f}) * ({dtheta_dt_val:.4f})")
    
    # Print the final result
    print(f"\nThe calculated vertical velocity of the bridge edge is: {vertical_velocity:.2f} m/min")

# Execute the function to solve the problem
solve_drawbridge_velocity()

# Final answer for the platform
final_answer = 37.79
# print(f'<<<{final_answer}>>>')