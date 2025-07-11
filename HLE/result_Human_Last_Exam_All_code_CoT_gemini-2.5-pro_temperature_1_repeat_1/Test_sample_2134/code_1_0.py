import math
import numpy as np

def solve_trajectory():
    """
    This function calculates the final position of the particle based on the problem description.
    """
    # Time at which to find the position
    t = 2 * math.sqrt(3)

    # The initial position x(0) is given by 3 + z, where z is a specific constant.
    # z = (18 - 6*sqrt(3))^(1/3) + (18 + 6*sqrt(3))^(1/3)
    # This constant z is the real root of the cubic equation z^3 - 18z - 36 = 0.
    # We can find the roots of this polynomial using numpy.
    cubic_coeffs = [1, 0, -18, -36]
    roots = np.roots(cubic_coeffs)

    # The real root corresponds to the value of z.
    real_roots = roots[np.isreal(roots)]
    z = real_roots[0].real

    # The initial position x(0)
    x0 = 3 + z

    # The physical model suggests a trajectory x(t) = x(0) - t^2 / 4.
    # We calculate the displacement term t^2 / 4.
    displacement = t**2 / 4

    # The final position x(t)
    xt = x0 - displacement

    # --- Output the results ---
    print("Based on the physics of the system, the particle's trajectory is modeled by:")
    print("x(t) = x(0) - t^2 / 4\n")

    print("Step 1: Evaluate the initial position x(0).")
    print(f"The constant z is the real root of z^3 - 18z - 36 = 0, which is z ≈ {z:.6f}")
    print(f"The initial position is x(0) = 3 + z = 3 + {z:.6f} = {x0:.6f}\n")

    print("Step 2: Calculate the displacement at t = 2*sqrt(3).")
    print(f"The displacement term is t^2 / 4 = (2*sqrt(3))^2 / 4 = 12 / 4 = {displacement:.1f}\n")

    print("Step 3: Calculate the final position x(t).")
    print("The final equation is: x(2*sqrt(3)) = x(0) - 3")
    print("Plugging in the numbers:")
    # We print the final equation with the calculated values
    print(f"{xt:.6f} = {x0:.6f} - {displacement:.1f}")

solve_trajectory()