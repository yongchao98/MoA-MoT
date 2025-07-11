import math

def calculate_toric_code_degeneracy(n, m):
    """
    Calculates the ground space degeneracy (GSD) of the toric code on a plane
    with n smooth holes and m rough holes.

    The degeneracy is given by 2^k, where k is the number of logical qubits.
    For a toric code on an infinite plane (topologically a disk, g=0, 1 boundary)
    with n smooth and m rough holes punched in it, the total number of boundaries
    is b = 1 + n + m.

    The number of logical qubits is k = 2g + b - 1 = 2*0 + (1 + n + m) - 1 = n + m.
    Therefore, GSD = 2^(n+m).
    """
    if n < 0 or m < 0:
        print("Number of holes cannot be negative.")
        return

    # Calculate the exponent, which is the number of logical qubits k
    exponent_k = n + m

    # Calculate the ground space degeneracy
    degeneracy = 2**exponent_k

    # Print the explanation and the final equation
    print(f"Given n = {n} smooth holes and m = {m} rough holes on a plane.")
    print(f"The number of logical qubits k = n + m.")
    print(f"So, k = {n} + {m} = {exponent_k}.")
    print(f"The ground space degeneracy is 2^k.")
    print(f"Resulting equation: 2^({exponent_k}) = {degeneracy}")

# Example usage with some values for n and m.
# You can change these values to test other cases.
num_smooth_holes = 3
num_rough_holes = 2

calculate_toric_code_degeneracy(num_smooth_holes, num_rough_holes)