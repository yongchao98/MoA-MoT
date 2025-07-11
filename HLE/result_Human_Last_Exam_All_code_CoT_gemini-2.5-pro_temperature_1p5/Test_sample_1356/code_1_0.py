def solve_and_display_period():
    """
    This function provides the step-by-step thinking process for solving the physics problem
    and presents the final derived formula for the period of motion in a clear format.
    """
    # These are the numerical constants derived from the formula T = 2 * pi * sqrt( (8*R) / (3*g) )
    multiplier = 2
    numerator_factor = 8
    denominator_factor = 3
    
    print("This script solves for the period of a pendulum-like system consisting of two disks on a rod.")
    print("Here is the plan used to derive the formula:")
    print("1. A coordinate system was defined with 'x' as the horizontal position of the top disk's center and 'theta' as the rod's angle from the vertical.")
    print("2. The Lagrangian (L = T - V) was formulated for the system.")
    print("   - T (Kinetic Energy) included translational and rotational energy for both disks.")
    print("   - V (Potential Energy) was based on the height of the second, hanging disk.")
    print("3. The two Euler-Lagrange equations of motion were derived from the Lagrangian.")
    print("4. These equations were linearized for small oscillations using small-angle approximations.")
    print("5. The resulting coupled linear differential equations were solved to find the equation for a simple harmonic oscillator, yielding the angular frequency (omega).")
    print("   - The derived squared angular frequency is: omega^2 = (3 * g) / (8 * R)")
    print("6. The period of motion 'T' is calculated as 2 * pi / omega.")
    
    print("\nThe final formula for the period T, showing each derived numerical component, is:")
    
    # We print the final formula in a human-readable string format.
    print(f"T = {multiplier} * pi * sqrt( ({numerator_factor} * R) / ({denominator_factor} * g) )")

# Execute the function to display the solution steps and the final formula.
solve_and_display_period()