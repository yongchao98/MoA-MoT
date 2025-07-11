import math

def calculate_gsd(n, m):
    """
    Calculates the ground space degeneracy (GSD) of the toric code with
    n smooth holes and m rough holes on a planar surface.

    The formula for the number of logical qubits, k, is k = n + m.
    The GSD is 2^k.

    Args:
        n (int): The number of smooth holes.
        m (int): The number of rough holes.
    """
    if not isinstance(n, int) or not isinstance(m, int) or n < 0 or m < 0:
        print("Error: The number of holes (n and m) must be non-negative integers.")
        return

    # For a toric code on a planar surface (a disk) with n smooth and m rough
    # internal holes, the number of logical qubits is k = n + m.
    num_logical_qubits = n + m

    # The ground space degeneracy is 2 raised to the power of the number of logical qubits.
    gsd = int(math.pow(2, num_logical_qubits))

    print(f"For a toric code with n={n} smooth holes and m={m} rough holes on a planar surface:")
    print("The number of logical qubits is k = n + m.")
    print(f"The ground space degeneracy is 2^k, which is calculated as:")
    # The final print statement shows each number in the equation as requested.
    print(f"2^({n} + {m}) = {gsd}")


# Example usage with n=2 smooth holes and m=3 rough holes.
# You can change these values to see other results.
n_smooth_holes = 2
m_rough_holes = 3

calculate_gsd(n_smooth_holes, m_rough_holes)