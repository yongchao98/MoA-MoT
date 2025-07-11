import math

def calculate_toric_code_degeneracy(m, n):
    """
    Calculates the ground space degeneracy (GSD) of the toric code
    with n smooth holes and m rough holes on a planar surface.

    The formula GSD = 2^(m + n - 2) is used, which is valid for m>=1 and n>=1.

    Args:
        m (int): The number of rough holes (must be >= 1).
        n (int): The number of smooth holes (must be >= 1).
    """
    if m < 1 or n < 1:
        print("This formula is intended for cases where m>=1 and n>=1.")
        print(f"The provided values are m={m}, n={n}.")
        # Handle specific edge cases for a more complete answer.
        if m + n == 0:
            degeneracy = 1 # Sphere
            print("GSD for a sphere (m=0, n=0) is 1.")
        elif m >= 1 and n == 0:
            degeneracy = 2**(m - 1) if m > 0 else 1
            print(f"GSD for m={m}, n=0 is {degeneracy}.")
        elif n >= 1 and m == 0:
            degeneracy = 2**(n - 1) if n > 0 else 1
            print(f"GSD for m=0, n={n} is {degeneracy}.")
        return

    # Calculate the exponent for the general formula
    exponent = m + n - 2
    
    # Calculate the ground space degeneracy
    degeneracy = 2**exponent
    
    # Print the explanation and the final equation with all numbers
    print(f"Calculating the Ground Space Degeneracy (GSD) for:")
    print(f"m = {m} (rough holes)")
    print(f"n = {n} (smooth holes)")
    print("\nUsing the formula for the general case (m>=1, n>=1):")
    # The prompt requires outputting each number in the final equation.
    print(f"GSD = 2^({m} + {n} - 2)")
    print(f"GSD = 2^{exponent}")
    print(f"GSD = {degeneracy}")

# --- Example Execution ---
# You can change these values to see the result for different numbers of holes.
num_rough_holes = 3
num_smooth_holes = 4

calculate_toric_code_degeneracy(num_rough_holes, num_smooth_holes)
