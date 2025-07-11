def calculate_toric_code_gsd(n, m):
    """
    Calculates the ground space degeneracy (GSD) of the toric code
    on a plane with n smooth and m rough holes.

    Args:
        n (int): The number of smooth holes.
        m (int): The number of rough holes.
    """
    # The number of logical qubits (k) is determined by the topology.
    # For a surface with genus g and b boundaries, k = 2*g + b - 1.
    # A plane has genus g=0.
    # The total number of boundaries b is n + m + 1 (including the boundary at infinity).
    # So, k = 2*0 + (n + m + 1) - 1 = n + m.
    
    k = n + m
    
    # The Ground Space Degeneracy (GSD) is 2^k.
    gsd = 2**k
    
    # Print the final equation with all numbers, as requested.
    # For example, for n=2 and m=3, this will print "2^(2+3) = 32"
    print(f"2^({n}+{m}) = {gsd}")

# Example usage with n=2 smooth holes and m=3 rough holes.
num_smooth_holes = 2
num_rough_holes = 3
calculate_toric_code_gsd(num_smooth_holes, num_rough_holes)