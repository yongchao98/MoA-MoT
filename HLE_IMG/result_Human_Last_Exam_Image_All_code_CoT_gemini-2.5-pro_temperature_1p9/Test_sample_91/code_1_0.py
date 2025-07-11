import math

def solve_drawbridge_velocity():
    """
    Calculates the vertical velocity of a drawbridge tip based on given parameters.
    """
    # Given parameters from the problem
    L = 40.0  # Length of the drawbridge in meters
    y_target = 10.0  # Target vertical height in meters
    
    # The given rate of change of theta is a constant value.
    # d(theta)/dt = - (3 * pi / 10) / cos(pi / 12)
    pi = math.pi
    cos_pi_over_12 = math.cos(pi / 12)
    dtheta_dt = -(3 * pi / 10) / cos_pi_over_12
    
    # Step 1: Explain the relationship between height y and angle theta
    print("The relationship between the vertical height (y), length (L), and angle (theta) is:")
    print("y = L * cos(theta)")
    print("\nThe formula for vertical velocity (v_y) is found by differentiating with respect to time (t):")
    print("v_y = dy/dt = -L * sin(theta) * d(theta)/dt")
    print("-" * 30)

    # Step 2: Find sin(theta) when y = 10m
    print(f"When the height y = {y_target} m and L = {L} m:")
    # Calculate cos(theta)
    cos_theta = y_target / L
    print(f"cos(theta) = y / L = {y_target} / {L} = {cos_theta}")
    
    # Calculate sin(theta)
    sin_theta = math.sqrt(1 - cos_theta**2)
    print(f"sin(theta) = sqrt(1 - cos(theta)^2) = sqrt(1 - {cos_theta**2:.4f}) = {sin_theta:.5f}")
    print("-" * 30)
    
    # Step 3: Present the given rate d(theta)/dt
    print("The rate of change of theta is given as a constant value:")
    print("d(theta)/dt = -(3 * pi) / (10 * cos(pi/12))")
    print(f"Numerically, this is: -(3 * {pi:.5f}) / (10 * {cos_pi_over_12:.5f}) = {dtheta_dt:.5f} rad/s")
    print("-" * 30)
    
    # Step 4: Assemble the final equation with all numbers
    print("Plugging all the values into the velocity equation:")
    print("v_y = -L * sin(theta) * d(theta)/dt")
    print(f"v_y = -({L}) * ({sin_theta:.5f}) * ({dtheta_dt:.5f})")
    
    # Perform the final calculation
    v_y = -L * sin_theta * dtheta_dt
    
    print("\nFinal calculated vertical velocity:")
    print(f"v_y = {v_y:.2f} m/s")

# Run the solver
solve_drawbridge_velocity()