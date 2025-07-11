def solve_icosahedron_puzzle():
    """
    This program explains the reasoning to determine the shape of the water
    surface in a half-filled icosahedron tank standing on one face.
    """
    print("Problem: What is the shape of the water surface in a half-filled icosahedron tank standing on a face?")
    print("-" * 80)
    print("Step 1: Analyzing the orientation of the icosahedron tank.")
    print("The tank stands on a face, meaning its base is a horizontal equilateral triangle.")
    print("The vertical axis passing through the center of this face is an axis of 3-fold rotational symmetry.\n")

    print("Step 2: Understanding the 'half-filled' condition.")
    print("A half-filled tank contains water occupying exactly half its total volume.")
    print("Due to the perfect central symmetry of the icosahedron, the horizontal plane that divides the volume in half must pass through the geometric center of the tank.\n")

    print("Step 3: Determining the shape of the water surface.")
    print("The water surface is the shape of the central, horizontal cross-section of the icosahedron.\n")

    print("Step 4: Using symmetry to deduce the final shape.")
    symmetries = {
        'rotational': 3,
        'central': 2
    }
    print(f" - The cross-section is perpendicular to an axis of {symmetries['rotational']}-fold symmetry, so the shape must have {symmetries['rotational']}-fold rotational symmetry.")
    print(f" - The cross-section passes through the center of the icosahedron, so it must also have point symmetry (or {symmetries['central']}-fold rotational symmetry).")
    
    final_symmetry = 6 # Least common multiple of 2 and 3
    print(f"\nEquation of Symmetry: A shape with both {symmetries['rotational']}-fold and {symmetries['central']}-fold symmetry must have {final_symmetry}-fold rotational symmetry.")

    print("\nStep 5: Conclusion.")
    print(f"A polygon that has {final_symmetry}-fold rotational symmetry is a regular hexagon.")
    print("\nTherefore, the shape of the water surface is a regular hexagon.")

solve_icosahedron_puzzle()
<<<A regular hexagon>>>