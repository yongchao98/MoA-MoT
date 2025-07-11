import math

def solve_for_epicycle_parameters():
    """
    Calculates the parameters (R, phi) for a deferent-epicycle model
    approximating motion around a square.
    """
    # Step 1: Determine the frequency ratio `phi`.
    # Based on the 4-fold rotational symmetry of the square path, the simplest
    # deferent-epicycle model that reproduces this symmetry requires a retrograde
    # epicycle with a frequency of -3 times the deferent's frequency.
    phi = -3

    # Step 2: Formulate the equation for the radius ratio `R`.
    # We match the ratio of max to min angular velocity of the model to that of the
    # object on the square, which is 2.
    # For a model with phi = -3, the angular velocity ratio is:
    # (R+3)(R+1) / ((R-1)(R-3)) = 2
    # This simplifies to the quadratic equation: R^2 - 12R + 3 = 0.
    
    # The coefficients of the quadratic equation a*R^2 + b*R + c = 0
    a = 1
    b = -12
    c = 3
    
    print("The problem reduces to solving the following quadratic equation for the radius ratio R:")
    print(f"{a}*R^2 + ({b})*R + {c} = 0")

    # Step 3: Solve the quadratic equation for R.
    # Using the quadratic formula: R = (-b +/- sqrt(b^2 - 4ac)) / (2a)
    discriminant = b**2 - 4 * a * c
    
    if discriminant < 0:
        print("The equation has no real solutions.")
        return

    root1 = (-b + math.sqrt(discriminant)) / (2 * a)
    root2 = (-b - math.sqrt(discriminant)) / (2 * a)

    # Step 4: Select the physically correct root.
    # For the object's angular velocity to always be positive, R must be > 3.
    if root1 > 3:
        R = root1
    elif root2 > 3:
        R = root2
    else:
        # This case should not be reached with the current equation.
        print("No physical solution found for R > 3.")
        return

    # Step 5: Output the final ordered pair (R, phi).
    print("\n--- Model Parameters ---")
    print(f"Radius Ratio (R = R_deferent / R_epicycle): {R}")
    print(f"Frequency Ratio (phi = omega_epicycle / omega_deferent): {phi}")
    print("\nThe ordered pair (R, phi) is:")
    print(f"({R}, {phi})")

solve_for_epicycle_parameters()