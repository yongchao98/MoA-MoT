import math

def calculate_toric_code_gsd(n, m):
    """
    Calculates the ground space degeneracy (GSD) of the toric code on a 
    2D plane with n smooth holes and m rough holes.

    Args:
        n (int): The number of smooth holes.
        m (int): The number of rough holes.
    """
    # Step 1: Explain the physical model.
    # The problem describes a toric code with n smooth and m rough holes.
    # A common and practical model for this is the "planar code", which resides on a 2D plane.
    # Topologically, a 2D plane is a sphere (genus g=0) with one hole (the boundary at infinity).
    # The n smooth and m rough holes are additional boundaries on this surface.
    print(f"Model: Toric code on a 2D plane with {n} smooth and {m} rough holes.")
    print("This corresponds to a topological surface with genus g=0.")

    # Step 2: Calculate the total number of boundaries (b).
    # Total boundaries = n (smooth) + m (rough) + 1 (at infinity).
    b = n + m + 1
    print(f"The total number of boundaries is b = n + m + 1 = {n} + {m} + 1 = {b}.")

    # Step 3: Calculate the number of logical qubits (k).
    # The formula for the number of logical qubits is k = 2*g + b - 1 (for b > 0).
    g = 0
    k = 2 * g + b - 1
    print(f"The number of logical qubits is k = 2*g + b - 1 = 2*({g}) + {b} - 1 = {k}.")
    
    # Step 4: Calculate the Ground Space Degeneracy (GSD).
    # GSD = 2^k
    gsd = int(math.pow(2, k))
    
    # Step 5: Output the final equation and result, showing each number.
    # The derived formula is GSD = 2^(n+m).
    print("\nThe final formula for the Ground Space Degeneracy (GSD) is 2^(n+m).")
    print("For the given values:")
    print(f"GSD = 2^({n} + {m}) = 2^{k} = {gsd}")

# Example values for the number of holes
n_example = 3  # Example number of smooth holes
m_example = 2  # Example number of rough holes

calculate_toric_code_gsd(n_example, m_example)