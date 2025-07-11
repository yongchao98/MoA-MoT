def calculate_toric_code_degeneracy(n, m):
    """
    Calculates the ground space degeneracy of the toric code with n smooth and m rough holes.

    Args:
        n (int): The number of smooth holes.
        m (int): The number of rough holes.
    """
    # The number of logical qubits 'k' is n + m based on the topological interpretation.
    k = n + m
    
    # The ground space degeneracy (GSD) is 2^k.
    degeneracy = 2**k
    
    # Print the explanation and the final equation with all numbers.
    print(f"For a toric code with n={n} smooth holes and m={m} rough holes:")
    print(f"The number of logical qubits is k = n + m = {n} + {m} = {k}.")
    print(f"The ground space degeneracy is GSD = 2^k.")
    print(f"GSD = 2^({n} + {m}) = {degeneracy}")


# Example values for demonstration
num_smooth_holes = 3
num_rough_holes = 2

calculate_toric_code_degeneracy(num_smooth_holes, num_rough_holes)