import math
from sympy.physics.wigner import wigner_3j
from sympy import N

def calculate_infinity_norm(j1, j2, j3):
    """
    Calculates the infinity-norm (maximum absolute value) for a Wigner 3j symbol
    by iterating through all valid m values.
    """
    max_abs_val = 0.0
    # Iterate over all possible m1 and m2 values
    for m1 in range(-j1, j1 + 1):
        for m2 in range(-j2, j2 + 1):
            m3 = -m1 - m2
            # The Wigner 3j symbol is non-zero only if |m3| <= j3.
            # The wigner_3j function handles this check automatically by returning 0.
            if abs(m3) <= j3:
                # Calculate the Wigner 3-j symbol and convert it to a float
                val = N(wigner_3j(j1, j2, j3, m1, m2, m3))
                abs_val = abs(float(val))
                if abs_val > max_abs_val:
                    max_abs_val = abs_val
    return max_abs_val

def solve_ratio():
    """
    Calculates the ratio of the maximum to the minimum infinity-norm for the
    nine Wigner 3-j symbols shown in the image.
    """
    # Define the (j1, j2, j3) parameters for each of the 9 plots.
    # j3 is fixed at 5. j1 and j2 range from 4 to 6.
    symbol_params = [
        (4, 4, 5), (4, 5, 5), (4, 6, 5),  # Row 1
        (5, 4, 5), (5, 5, 5), (5, 6, 5),  # Row 2
        (6, 4, 5), (6, 5, 5), (6, 6, 5)   # Row 3
    ]

    # Calculate the infinity-norm for each symbol
    infinity_norms = []
    for params in symbol_params:
        j1, j2, j3 = params
        norm = calculate_infinity_norm(j1, j2, j3)
        infinity_norms.append(norm)

    # Find the maximum and minimum of these norms
    # Filter out any norms that are zero to avoid division by zero
    valid_norms = [norm for norm in infinity_norms if norm > 0]
    if not valid_norms:
        print("No non-zero norms found.")
        return

    max_norm = max(valid_norms)
    min_norm = min(valid_norms)

    # Calculate the ratio
    ratio = max_norm / min_norm

    # Print the final results in the required format
    print(f"Maximum infinity-norm: {max_norm}")
    print(f"Minimum infinity-norm: {min_norm}")
    print(f"Ratio = {max_norm} / {min_norm} = {ratio}")

if __name__ == '__main__':
    solve_ratio()
