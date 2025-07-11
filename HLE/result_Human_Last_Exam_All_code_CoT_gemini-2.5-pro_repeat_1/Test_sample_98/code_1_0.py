import math

def solve_icosahedron_puzzle():
    """
    Deduces the shape of the water surface in a half-filled icosahedron
    tank standing on one of its faces.
    """

    # --- Properties of the Icosahedron and the Setup ---
    
    # 1. An icosahedron is a regular polyhedron with 20 faces.
    num_faces = 20
    
    # 2. It's "half-filled", meaning the volume of water is 1/2 the total volume.
    fill_fraction_numerator = 1
    fill_fraction_denominator = 2
    
    # 3. When an icosahedron stands on a face, the axis perpendicular to that
    #    face is an axis of 3-fold rotational symmetry.
    rotational_symmetry_fold = 3

    # --- Logical Deduction ---

    print("Step 1: Locating the water surface")
    print("A regular icosahedron has perfect point symmetry about its center.")
    print(f"Any plane passing through the center divides the volume into two equal halves ({fill_fraction_numerator}/{fill_fraction_denominator}).")
    print("Since the tank is half-filled, the water surface must be a horizontal plane passing through the icosahedron's center.")
    print("-" * 50)

    print("Step 2: Determining the shape from symmetry")
    print("The shape of the surface is the shape of the central cross-section of the icosahedron.")
    print(f"The axis perpendicular to the base is an axis of {rotational_symmetry_fold}-fold rotational symmetry.")
    print(f"This means the cross-section shape must look the same if rotated by 360 / {rotational_symmetry_fold} = {360 // rotational_symmetry_fold} degrees.")
    
    print("\nAdditionally, the icosahedron has central (or inversion) symmetry.")
    print("Because the cutting plane passes through the center, the 2D shape of the cross-section must also be centrally symmetric.")
    print("-" * 50)

    print("Conclusion: The Final Shape")
    print(f"A 2D shape that has both {rotational_symmetry_fold}-fold rotational symmetry and central symmetry must be a regular hexagon.")
    print("A hexagon has 6 sides, and being regular means all sides and angles are equal.")
    
    final_shape = "regular hexagon"
    
    print(f"\nTherefore, the shape of the water surface is a {final_shape}.")

# Execute the puzzle solver
solve_icosahedron_puzzle()