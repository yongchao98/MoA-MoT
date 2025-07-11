import math

def solve_projectile_problem():
    """
    Calculates the initial speed 'u' for the rock and the memory usage 'z'.
    """
    # Step 1: Define the problem parameters
    initial_distance = 300.0  # meters
    lion_speed = 5.0      # m/s
    launch_angle_deg = 60.0 # degrees
    gravity = 9.8         # m/s^2

    # Step 2: Derive the quadratic equation for the initial speed 'u'.
    # The physics of the problem leads to a quadratic equation for u:
    # (sin(2*a)) * u^2 + (2 * v * sin(a)) * u - (d * g) = 0
    # where a is the launch angle, v is the lion's speed, d is the initial distance.
    # This equation is in the form A*u^2 + B*u - C = 0.

    # Step 3: Calculate the coefficients A, B, and C for the equation.
    launch_angle_rad = math.radians(launch_angle_deg)
    sin_a = math.sin(launch_angle_rad)
    sin_2a = math.sin(2 * launch_angle_rad)

    A = sin_2a
    B = 2 * lion_speed * sin_a
    C = initial_distance * gravity

    print("The problem is solved by finding the positive root of the quadratic equation: A*u^2 + B*u - C = 0")
    print("\nBased on the given values:")
    print(f"  - Initial distance (d) = {initial_distance} m")
    print(f"  - Lion's speed (v) = {lion_speed} m/s")
    print(f"  - Launch angle (a) = {launch_angle_deg} degrees")
    print(f"  - Gravity (g) = {gravity} m/s^2")

    print("\nThe coefficients of the equation are:")
    print(f"A = sin(2 * {launch_angle_deg}°) = {A:.4f}")
    print(f"B = 2 * {lion_speed} * sin({launch_angle_deg}°) = {B:.4f}")
    print(f"C = {initial_distance} * {gravity} = {C:.1f}")

    # Step 4: Display the final equation and solve for u.
    # We solve for u using the quadratic formula: u = [-B + sqrt(B^2 + 4AC)] / 2A
    print("\nThe final equation to solve for u is:")
    print(f"({A:.4f})*u^2 + ({B:.4f})*u - {C:.1f} = 0")
    
    discriminant = B**2 + 4 * A * C
    u_speed = (-B + math.sqrt(discriminant)) / (2 * A)
    u_rounded = round(u_speed, 1)

    print(f"\nThe required initial speed 'u' is {u_speed:.4f} m/s, which rounds to {u_rounded} m/s.")

    # Step 5: Calculate the memory usage 'z'.
    # A hypothetical program would need variables for the core inputs (d, v)
    # and the result (u). We assume these would be of the 'frac' type to handle
    # non-integer math.
    # Size of 'frac' = sizeof(n: signed char) + sizeof(d: unsigned char) + sizeof(e: signed char)
    # Size of 'char' is 2D.
    # Size of 'frac' = 2D + 2D + 2D = 6D.
    # Number of variables = 3 (for d, v, u).
    # Memory usage 'z' = 3 * 6D = 18D.
    z_memory = 18
    print(f"The memory usage 'z' for the program's variables is {z_memory} D.")
    
    # Step 6: Print the final answer in the required format.
    print("\nFinal Answer (u:z):")
    print(f"<<<{u_rounded}:{z_memory}>>>")

solve_projectile_problem()