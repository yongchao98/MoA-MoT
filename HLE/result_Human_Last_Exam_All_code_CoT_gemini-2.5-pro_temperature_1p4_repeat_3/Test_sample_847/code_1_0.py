import numpy as np
import math

def solve_and_find_cost():
    """
    This function solves for the optimal strategy in the sorting game
    and calculates the minimal cost coefficient.
    """
    # Step 1: Define the cubic equation p^3 + p - 1 = 0.
    # The coefficients are for p^3, p^2, p^1, p^0.
    coeffs = [1, 0, 1, -1]

    # Step 2: Solve the equation. np.roots finds all roots (real and complex).
    roots = np.roots(coeffs)

    # We need the single real root, which corresponds to a valid probability.
    p0 = roots[np.isreal(roots)].real[0]

    # Step 3: Use the simplified formula to calculate the minimal cost per bit (C_min).
    # C_min = -1 / log2(p0)
    log2_p0 = math.log2(p0)
    C_min = -1 / log2_p0

    # Step 4: Print the results, showing the numbers in the final equation as requested.
    print("The optimal strategy involves asking questions where the 'yes' probability, p0, is the real root of the equation p^3 + p - 1 = 0.")
    print(f"The value of this probability is p0 ≈ {p0:.4f}")
    print("\nThe minimal cost per bit of information, C_min, is given by the formula:")
    print("C_min = -1 / log2(p0)")
    print("\nPlugging in the numbers:")
    print(f"log2({p0:.4f}) ≈ {log2_p0:.4f}")
    print(f"C_min = -1 / ({log2_p0:.4f}) ≈ {C_min:.4f}")

    # Step 5: Print the final answer, rounded to 3 decimal places.
    print(f"\nThe minimal cost coefficient, rounded to 3 decimal places, is {C_min:.3f}")

solve_and_find_cost()