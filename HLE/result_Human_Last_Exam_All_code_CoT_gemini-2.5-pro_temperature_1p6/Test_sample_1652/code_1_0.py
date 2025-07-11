import math

def solve_for_u():
    """
    Calculates the initial speed 'u' for the rock to hit the lion and
    the memory usage 'z' on the Wuxing architecture.
    """

    # --- Part 1: Physics Calculation ---
    # Initial conditions from the problem statement
    d_lion = 300  # meters
    v = 5         # m/s
    a_deg = 60    # degrees
    g = 9.8       # m/s^2

    # The physics of the problem can be modeled with a quadratic equation of the form:
    # A*u^2 + B*u + C = 0
    # This is derived from setting the rock's range equal to the lion's final position.
    # sin(2*a)*u^2 + (2*v*sin(a))*u - d*g = 0
    
    # Convert angle to radians for Python's math functions
    a_rad = math.radians(a_deg)

    # Calculate the coefficients of the quadratic equation
    coeff_A = math.sin(2 * a_rad)
    coeff_B = 2 * v * math.sin(a_rad)
    coeff_C = -d_lion * g

    # Solve for u using the quadratic formula: u = (-B ± sqrt(B^2 - 4AC)) / 2A
    # We take the positive root, as speed cannot be negative.
    discriminant = (coeff_B**2) - 4 * coeff_A * coeff_C
    
    if discriminant < 0:
        print("No real solution exists for the initial speed u.")
        return

    u = (-coeff_B + math.sqrt(discriminant)) / (2 * coeff_A)
    
    # Round the result for 'u' to one decimal place
    u_rounded = round(u, 1)

    # --- Part 2: Memory Usage Calculation (z) ---
    # As per the Wuxing architecture:
    # size of int = 5D
    # size of frac = size(n) + size(d) + size(e) = 2D + 2D + 2D = 6D

    # We determine the data types for the variables in the program:
    # d_lion (300): int
    # v (5): int
    # a (60): int
    # rock_mass (0.5 kg, from problem): frac
    # g (9.8): frac
    # u (the result, a non-integer): frac
    
    mem_int_count = 3
    mem_frac_count = 3
    
    size_int = 5
    size_frac = 6
    
    z = (mem_int_count * size_int) + (mem_frac_count * size_frac)

    # --- Final Output ---
    # As requested, outputting each number in the final equation.
    # The equation is A*u^2 + B*u + C = 0
    print("The final equation is derived from the projectile motion formulas.")
    print("It takes the form A*u^2 + B*u + C = 0, where:")
    print(f"A (from sin(2a)): {coeff_A}")
    print(f"B (from 2*v*sin(a)): {coeff_B}")
    print(f"C (from -d*g): {coeff_C}")
    print("-" * 20)
    
    # Print the final result in the format u:z
    print(f"Result in the format u:z is {u_rounded}:{z}")

solve_for_u()