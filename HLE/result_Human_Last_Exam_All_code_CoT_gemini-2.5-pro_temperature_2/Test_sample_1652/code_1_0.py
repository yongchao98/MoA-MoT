import math

def solve_problem():
    """
    Solves the physics problem and calculates memory usage.
    """
    # Part 1: Physics Calculation for u
    
    # Given values
    d = 300.0  # initial distance in meters
    v = 5.0    # lion's speed in m/s
    angle_deg = 60.0
    g = 9.8    # acceleration due to gravity in m/s^2

    # Convert angle to radians for math functions
    a_rad = math.radians(angle_deg)

    # The problem can be modeled with the following quadratic equation for u:
    # (sin(2a)/g) * u^2 + (2*v*sin(a)/g) * u - d = 0
    # Let's define the coefficients A, B, C for Au^2 + Bu + C = 0
    
    A = math.sin(2 * a_rad)
    B = 2 * v * math.sin(a_rad)
    C = -d * g

    # Using the quadratic formula: u = [-B ± sqrt(B^2 - 4AC)] / 2A
    # We use the adjusted coefficients where we multiply the whole equation by g
    # A_u^2 + B_u + C = 0
    
    discriminant = B**2 - 4 * A * C
    
    # Speed u must be positive, so we take the positive root
    u = (-B + math.sqrt(discriminant)) / (2 * A)
    
    u_rounded = round(u, 1)

    # Part 2: Memory Calculation for z

    # Based on the Wuxing architecture specification:
    # int = 5D
    # frac = 6D (signed char (2D) + unsigned char (2D) + signed char (2D))

    # Variables from the problem description:
    # 1. Rock mass (0.5 kg): Non-integer, requires 'frac'.
    size_mass = 6  # frac mass = 1/2;

    # 2. Initial distance (300 m): Integer, fits in 'int'.
    size_distance = 5  # int d = 300;

    # 3. Lion speed (5 m/s): Integer, fits in 'int'.
    size_lion_speed = 5  # int v = 5;

    # 4. Angle (60 degrees): Integer, fits in 'int'.
    size_angle = 5  # int a = 60;

    # 5. Gravity (9.8 m/s^2): Non-integer, requires 'frac'.
    size_gravity = 6  # frac g = 98/10;
    
    # 6. Initial speed u (result): Will be a non-integer, requires 'frac'.
    size_u = 6 # frac u;

    # Total memory usage z
    z = size_mass + size_distance + size_lion_speed + size_angle + size_gravity + size_u
    
    print(f"{u_rounded}:{z}")

solve_problem()