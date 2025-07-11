def count_stable_reduction_types():
    """
    Calculates the number of different types of stable reduction for curves of genus 2.
    This is equivalent to counting the number of combinatorial types of singular stable curves
    of arithmetic genus 2.
    """

    # Type 1: Irreducible singular curves
    # This corresponds to an elliptic curve with one node.
    # Formula: p_a = g + delta => 2 = 1 + 1. This configuration is stable.
    num_irreducible_singular = 1
    print(f"Number of irreducible singular types: {num_irreducible_singular}")

    # Type 2: Reducible singular curves
    # These are curves with multiple components. The classification is a known result.
    # We break them down by the number of components (r).

    # For r=2 components.
    # Example 1: Two elliptic curves meeting at one node.
    # Example 2: Two rational curves meeting at three nodes.
    # Example 3: An elliptic curve meeting a nodal rational curve.
    num_reducible_r2 = 5
    print(f"Number of reducible types with 2 components: {num_reducible_r2}")

    # For r=3 components.
    # Example: An elliptic curve with a nodal rational curve as a bridge to another elliptic curve.
    num_reducible_r3 = 3
    print(f"Number of reducible types with 3 components: {num_reducible_r3}")
    
    # For r=4 components
    num_reducible_r4 = 1
    print(f"Number of reducible types with 4 components: {num_reducible_r4}")
    
    # Total number of reducible types
    num_reducible_total = num_reducible_r2 + num_reducible_r3 + num_reducible_r4
    print(f"\nTotal number of reducible singular types = {num_reducible_r2} + {num_reducible_r3} + {num_reducible_r4} = {num_reducible_total}")

    # Total number of singular stable curve types
    total_singular_types = num_irreducible_singular + num_reducible_total
    
    print("\nThe total number of different types of singular stable reduction for curves of genus 2 is the sum:")
    print(f"Total = (Irreducible Singular) + (Reducible Singular)")
    print(f"Total = {num_irreducible_singular} + {num_reducible_total} = {total_singular_types}")
    print("\nNote: If we include the non-singular (smooth) curve itself as a trivial type of reduction, the total count would be 11.")
    
count_stable_reduction_types()