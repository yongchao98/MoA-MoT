def check_homotopy_section_for_surface(genus, num_boundaries, is_orientable=True):
    """
    Checks if the fibration pi_{k,l} admits a homotopy section for the interior of a
    compact surface with the given properties.

    The existence condition is as follows:
    1. If the manifold is non-compact (i.e., the original surface had boundaries),
       a section always exists.
    2. If the manifold is compact and closed (i.e., the original surface had no boundaries),
       a homotopy section exists if and only if its Euler characteristic is zero.

    Args:
        genus (int): The genus of the surface.
        num_boundaries (int): The number of boundary components of the compact surface.
        is_orientable (bool): Whether the surface is orientable.
    """
    # The manifold M is the interior of a compact surface.
    # If the surface has boundaries, its interior M is non-compact.
    if num_boundaries > 0:
        print(f"The manifold is non-compact because the number of boundaries is {num_boundaries}, which is greater than 0.")
        print("A homotopy section exists.")
        print(f"\nFinal check: number of boundaries = {num_boundaries} > 0")

    # If the surface has no boundaries, M is a closed manifold.
    else:
        print(f"The manifold is closed (compact, no boundary) because the number of boundaries is {num_boundaries}.")
        if is_orientable:
            # For an orientable closed surface, the Euler characteristic is chi = 2 - 2g.
            euler_char = 2 - 2 * genus
            print(f"The Euler characteristic is calculated as: chi = 2 - 2 * g = 2 - 2 * {genus} = {euler_char}.")
            print(f"\nFinal equation: 2 - 2 * {genus} = {euler_char}")

        else:
            # For a non-orientable closed surface, the Euler characteristic is chi = 2 - g.
            euler_char = 2 - genus
            print(f"The Euler characteristic is calculated as: chi = 2 - g = 2 - {genus} = {euler_char}.")
            print(f"\nFinal equation: 2 - {genus} = {euler_char}")

        if euler_char == 0:
            print("A homotopy section exists because the Euler characteristic is 0.")
        else:
            print("A homotopy section does not exist because the Euler characteristic is not 0.")


# --- Example 1: The 2-Sphere (satisfies A and C, but no section exists) ---
print("--- Test Case 1: Sphere (M = S^2) ---")
# M is interior of a ball with its boundary collapsed to a point, i.e., a sphere.
# This is a closed manifold (genus=0, boundaries=0).
check_homotopy_section_for_surface(genus=0, num_boundaries=0, is_orientable=True)

# --- Example 2: The Torus (violates B, but a section exists) ---
print("\n--- Test Case 2: Torus (M = T^2) ---")
# This is a closed manifold (genus=1, boundaries=0).
check_homotopy_section_for_surface(genus=1, num_boundaries=0, is_orientable=True)

# --- Example 3: Interior of a Disk (satisfies B) ---
print("\n--- Test Case 3: Open Disk (M = R^2) ---")
# M is the interior of a compact disk (genus=0, boundaries=1).
check_homotopy_section_for_surface(genus=0, num_boundaries=1, is_orientable=True)