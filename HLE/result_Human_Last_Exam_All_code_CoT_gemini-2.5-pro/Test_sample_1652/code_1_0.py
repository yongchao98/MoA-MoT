import math

def solve_for_initial_speed():
    """
    Calculates the initial speed 'u' for the rock to hit the lion.
    """
    # --- Problem Parameters ---
    distance = 300.0  # Initial distance in meters
    lion_speed = 5.0  # Lion's speed in m/s
    gravity = 9.8     # Acceleration due to gravity in m/s^2
    angle_deg = 60.0  # Angle of throw in degrees

    # Convert angle to radians for use in Python's math functions
    angle_rad = math.radians(angle_deg)

    # --- Derive Quadratic Equation Coefficients ---
    # The physics of the problem leads to a quadratic equation for the initial speed u:
    # A*u^2 + B*u + C = 0
    # where:
    # A = sin(2*a)
    # B = 2 * v * sin(a)
    # C = -distance * g
    
    A = math.sin(2 * angle_rad)
    B = 2 * lion_speed * math.sin(angle_rad)
    C = -distance * gravity
    
    # --- Solve for u using the quadratic formula ---
    # u = (-B + sqrt(B^2 - 4AC)) / (2A)
    # We take the positive root because speed must be positive.
    discriminant = B**2 - 4 * A * C
    
    if discriminant < 0:
        print("The rock can never hit the lion under these conditions.")
        return

    u = (-B + math.sqrt(discriminant)) / (2 * A)

    # --- Print the results as requested ---
    print("To find the initial speed u, we solve the quadratic equation: A*u^2 + B*u + C = 0")
    print("\nBased on the given parameters:")
    print(f"d = {distance} m")
    print(f"v = {lion_speed} m/s")
    print(f"a = {angle_deg} degrees")
    
    print("\nThe coefficients of the equation are:")
    print(f"A = sin(2 * {angle_deg}°) = {A:.4f}")
    print(f"B = 2 * {lion_speed} * sin({angle_deg}°) = {B:.4f}")
    print(f"C = -{distance} * {gravity} = {C:.1f}")

    print("\nThe final equation with numerical values is:")
    print(f"({A:.4f})*u^2 + ({B:.4f})*u + ({C:.1f}) = 0")

    print("\nSolving for u:")
    print(f"u = (-({B:.4f}) + sqrt(({B:.4f})^2 - 4*({A:.4f})*({C:.1f}))) / (2*({A:.4f}))")
    
    print("\nThe required initial speed of the rock is:")
    print(f"u = {u:.1f} m/s")

solve_for_initial_speed()

<<<57.6:34>>>