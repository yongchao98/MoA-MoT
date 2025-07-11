import math

def solve_radiation_ratio():
    """
    This script calculates the maximum achievable ratio of bidirectional conical power
    to line intensity for a rotating charged sphere, as described in the problem.

    The derivation shows that the leading order radiation is from a rotating magnetic dipole,
    with an intensity pattern proportional to (1 + cos^2(theta)), where theta is the angle
    from the rotation axis.

    The maximum ratio is found by aligning the cone with the rotation axis and the line
    perpendicular to it. This leads to the analytical formula for the ratio:
    Ratio = (π/3) * (16 - 7√2)

    This script calculates the value of this expression.
    """

    # --- Step 1: Define the components from the final analytical equation ---
    # The equation is: Ratio = (pi / 3) * (16 - 7 * sqrt(2))
    pi_val = math.pi
    denominator = 3.0
    val_16 = 16.0
    val_7 = 7.0
    sqrt2_val = math.sqrt(2)

    # --- Step 2: Print each number in the final equation ---
    print("This script calculates the result from the derived formula: Ratio = (π / 3) * (16 - 7 * √2)\n")
    print("--- Individual Numbers in the Final Equation ---")
    print(f"Value for π (pi): {pi_val}")
    print(f"Value for the denominator '3': {denominator}")
    print(f"Value for the number '16': {val_16}")
    print(f"Value for the number '7': {val_7}")
    print(f"Value for √2 (sqrt(2)): {sqrt2_val}\n")

    # --- Step 3: Calculate the result step-by-step ---
    term_in_parentheses = val_16 - val_7 * sqrt2_val
    final_ratio = (pi_val / denominator) * term_in_parentheses

    # --- Step 4: Output the final answer ---
    print("--- Final Calculation ---")
    print(f"The calculated value of (16 - 7 * √2) is: {term_in_parentheses}")
    print(f"The maximum achievable ratio is (π / 3) * {term_in_parentheses} = {final_ratio}")

if __name__ == '__main__':
    solve_radiation_ratio()