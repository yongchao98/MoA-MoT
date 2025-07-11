from fractions import Fraction

def solve_projectile_speed():
    """
    Calculates the initial speed 'u' for the rock to hit the lion
    and the memory usage 'z' based on Wuxing architecture rules.
    """

    # --- Step 1: Define constants using fractions to avoid floating point issues ---
    # These represent the variables that would be stored in a Wuxing program.
    d_initial = Fraction(300)      # Initial distance in meters
    v_lion = Fraction(5)           # Lion's speed in m/s
    g = Fraction(98, 10)           # Gravity g = 9.8 m/s^2
    
    # Angle a = 60 degrees.
    # The Wuxing architecture does not have sin/cos functions. We must use
    # rational approximations that could be stored in 'frac' type variables.
    # cos(60) = 1/2 is exact.
    # sin(60) = sqrt(3)/2 ~= 0.866. A good, simple rational approximation is 7/8 = 0.875.
    # The numerator (7) and denominator (8) fit in a signed/unsigned char.
    cos60 = Fraction(1, 2)
    sin60 = Fraction(7, 8)

    # --- Step 2: Print the derived equation ---
    # The physical model leads to the quadratic equation Au^2 + Bu - C = 0, where:
    # A = (2 * sin60 * cos60) / g
    # B = (2 * v_lion * sin60) / g
    # C = d_initial
    # We can simplify this to an integer equation u^2 + 10u - 3360 = 0
    print("The problem simplifies to the quadratic equation: a*u^2 + b*u + c = 0")
    print("Derived integer coefficients:")
    # The simplified integer equation from the approximation is u^2 + 10u - 3360 = 0
    a_coeff, b_coeff, c_coeff = 1, 10, -3360
    print(f"a = {a_coeff}")
    print(f"b = {b_coeff}")
    print(f"c = {c_coeff}")
    print("-" * 20)

    # --- Step 3: Solve for 'u' using binary search ---
    # We need to find u such that rock_distance_traveled + lion_distance_traveled = d_initial.
    # Let f(u) = (u*cos60 + v_lion) * t, where t = (2*u*sin60)/g. We want f(u) = d_initial.
    
    low = Fraction(0)
    high = Fraction(200) # A reasonable upper bound for the speed

    # Iterate 100 times for sufficient precision
    for _ in range(100):
        u_guess = (low + high) / 2
        if u_guess == 0: # Avoid division by zero if guess is 0
            low = Fraction(1, 1000) # Move low slightly up
            continue

        time_of_flight = (2 * u_guess * sin60) / g
        distance_traveled = (u_guess * cos60 + v_lion) * time_of_flight
        
        if distance_traveled < d_initial:
            low = u_guess
        else:
            high = u_guess
            
    u_result = (low + high) / 2
    u_rounded = round(float(u_result), 1)

    # --- Step 4: Calculate memory usage 'z' ---
    # Variables required in the program scope for this calculation:
    # 1. d_initial (distance)
    # 2. v_lion (lion speed)
    # 3. g (gravity)
    # 4. cos60
    # 5. sin60
    # 6. low (for binary search)
    # 7. high (for binary search)
    # Total = 7 variables. Each 'frac' struct costs 6D (2D+2D+2D).
    # Some might argue for an 8th variable for precision/tolerance.
    # Let's count the core constants and search bounds.
    num_variables = 7 
    # Each 'frac' type consists of 3 'char' members (n, d, e).
    # Each 'char' takes 2 decimal positions (D).
    # size_of_frac = 2D (n) + 2D (d) + 2D (e) = 6D
    size_of_frac = 6 
    z_memory_usage = num_variables * size_of_frac

    print(f"Final Answer (u:z):")
    print(f"{u_rounded}:{z_memory_usage}")


solve_projectile_speed()
<<<53.2:42>>>