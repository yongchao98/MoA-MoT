import math

def solve_for_u():
    """
    Calculates the initial speed 'u' for the rock to hit the lion.
    """
    # 1. Define constants for the physics problem
    g = 9.8  # acceleration due to gravity in m/s^2
    v = 5.0    # speed of the lion in m/s
    dist = 300.0 # initial distance in meters
    angle_deg = 60.0 # launch angle in degrees

    # Convert angle to radians for math functions
    angle_rad = math.radians(angle_deg)

    # 2. Derive coefficients for the quadratic equation: A*u^2 + B*u + C = 0
    # From the physics equation: sin(2a)*u^2 + (2*v*sin(a))*u - dist*g = 0
    A = math.sin(2 * angle_rad)
    B = 2 * v * math.sin(angle_rad)
    C = -dist * g

    # As requested, print the equation with its numerical coefficients
    print("Derived equation to solve for u:")
    # The '+' before C is to correctly display its negative sign
    print(f"{A:.4f} * u^2 + {B:.4f} * u + ({C:.4f}) = 0")
    print("-" * 20)

    # 3. Solve the quadratic equation for u using the bisection method
    # This method is used because the Wuxing architecture lacks a sqrt function.
    def f(u):
        return A * u**2 + B * u + C

    u_low = 0.0
    u_high = 200.0  # An initial guess for the upper bound of the speed
    
    # Check if a solution exists in the interval
    if f(u_low) * f(u_high) >= 0:
        print("Bisection method cannot guarantee a solution in the initial interval.")
        return

    # Iterate to find the root
    for _ in range(100): # 100 iterations for high precision
        u_mid = (u_low + u_high) / 2
        f_mid = f(u_mid)

        if f_mid == 0:
            break
        elif f(u_low) * f_mid < 0:
            u_high = u_mid
        else:
            u_low = u_mid
            
    u_result = (u_low + u_high) / 2
    u_rounded = round(u_result, 1)

    # 4. Calculate memory usage 'z' for the Wuxing architecture
    # The size of one 'frac' variable is 2D(n) + 2D(d) + 2D(e) = 6D.
    # Variables needed for the solver: g, v, dist, sin(a), sin(2a), u_low, u_high
    num_variables = 7
    size_of_frac = 6 # in Decimal digits (D)
    z = num_variables * size_of_frac

    # 5. Print the final answer in the required format u:z
    print(f"Calculated result: u = {u_result:.4f} m/s")
    print(f"Memory usage: z = {z} D")
    print("-" * 20)
    print("Final answer:")
    print(f"{u_rounded}:{z}")


solve_for_u()