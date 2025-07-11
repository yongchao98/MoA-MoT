import math

def solve_for_speed():
    """
    Calculates the initial speed 'u' for the rock to hit the lion
    using the Bisection Method.
    """
    # 1. Define constants from the problem
    v = 5.0  # Lion's speed in m/s
    g = 9.8  # Gravitational acceleration in m/s^2
    distance = 300.0  # Initial distance in meters
    angle_deg = 60.0
    angle_rad = math.radians(angle_deg)
    sin_a = math.sin(angle_rad)

    # 2. Define the quadratic equation f(u) = u^2 + B*u + C = 0
    # derived from the physics.
    # B = 2*v
    # C = -distance*g / sin(a)
    B = 2 * v
    C = -distance * g / sin_a

    # The equation to solve is u^2 + 10*u - 3394.8... = 0
    
    # 3. Use Bisection Method to find the root 'u'
    low_u = 0.0
    high_u = 100.0  # A reasonable upper bound for the speed
    
    # The equation we are solving for u: f(u) = u^2 + B*u + C
    def f(u):
        return u**2 + B*u + C

    # Iterate to find the root with high precision
    for _ in range(100):
        mid_u = (low_u + high_u) / 2
        if f(mid_u) == 0:
            break
        # If f(mid_u) and f(low_u) have different signs, the root is in the lower half
        elif f(mid_u) * f(low_u) < 0:
            high_u = mid_u
        # Otherwise, the root is in the upper half
        else:
            low_u = mid_u
            
    # The final answer for u is the midpoint of the final range
    u = (low_u + high_u) / 2

    # As requested, output each number in the final equation being solved: u^2 + B*u + C = 0
    print(f"The equation being solved is: u^2 + ({B:.1f})*u + ({C:.1f}) = 0")
    print(f"When solved for u, the rock's initial speed must be: {u:.1f} m/s")

solve_for_speed()