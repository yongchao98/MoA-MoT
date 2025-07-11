import math
from fractions import Fraction

def solve_projectile_speed():
    """
    This function calculates the initial speed of the rock required to hit the lion,
    adhering to the constraints of the Wuxing computer architecture.
    """

    # Step 1: Define problem constants as fractions to simulate the 'frac' type.
    # g = 9.8 m/s^2 = 98/10 = 49/5
    g = Fraction(49, 5)
    # v = 5 m/s (lion's speed)
    v = Fraction(5, 1)
    # Initial distance = 300 m
    dist = Fraction(300, 1)
    # Angle a = 60 degrees. We need its sine and cosine.
    # cos(60) = 1/2 (exact)
    cos_a = Fraction(1, 2)
    # sin(60) = sqrt(3)/2. We use a rational approximation 13/15 (~0.8667)
    # as sqrt is not available.
    sin_a = Fraction(13, 15)

    # Step 2: Formulate the quadratic equation Au^2 + Bu + C = 0.
    # This is derived from the projectile motion equations:
    # t_flight = (2 * u * sin(a)) / g
    # dist = x_rock + x_lion = (u*cos(a)*t_flight) + (v*t_flight)
    # Substituting t_flight and rearranging gives:
    # (sin(a)cos(a)) * u^2 + (v*sin(a)) * u - (dist/2 * g) = 0
    A = sin_a * cos_a
    B = v * sin_a
    C = -(dist / 2) * g

    # Step 3: Solve the equation using the bisection method to avoid sqrt.
    # Define the function f(u) = Au^2 + Bu + C
    def f(u_val):
        return A * u_val**2 + B * u_val + C

    # Set initial bounds for the search. Speed u must be positive.
    low = Fraction(0)
    high = Fraction(200)  # A reasonable upper estimate for the speed

    # Iterate 100 times for high precision.
    for _ in range(100):
        mid = (low + high) / 2
        if f(mid) < 0:
            low = mid
        else:
            high = mid
    
    # The result is the final value of 'low' or 'high'.
    u_frac = low

    # Step 4: Calculate the final values for the answer.
    # u: the speed, rounded to 0.1
    u_rounded = round(float(u_frac), 1)

    # z: memory usage in D.
    # 6 variables (v, g, dist, sin_a, cos_a, u) are stored.
    # Each 'frac' type costs 6D (2D for n, 2D for d, 2D for e).
    num_variables = 6
    memory_per_variable = 6
    z_memory = num_variables * memory_per_variable

    # Step 5: Output the results as required.
    # First, output the numbers in the final equation being solved.
    print(f"Equation: ({A.numerator}/{A.denominator}) * u^2 + ({B.numerator}/{B.denominator}) * u + ({C.numerator}) = 0")
    
    # Then, print the final answer in the specified u:z format.
    print(f"Answer: {u_rounded}:{z_memory}")

solve_projectile_speed()