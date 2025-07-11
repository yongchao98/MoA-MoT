def solve_icosahedron_puzzle():
    """
    This function explains the reasoning to determine the shape of the water
    surface in a half-filled icosahedron tank standing on a face.
    """

    print("Analyzing the problem step by step:")
    print("------------------------------------")

    # Step 1: Define the setup
    print("Step 1: The tank is a regular icosahedron (20 triangular faces) standing on one face as its horizontal base.")

    # Step 2: Interpret "half-filled" using symmetry
    print("Step 2: A regular icosahedron has a center of symmetry. A plane passing through this center divides its volume into two equal halves.")
    print("        Therefore, for a half-filled tank, the water surface must be a horizontal plane passing through the icosahedron's center.")

    # Step 3: Identify the geometric shape of the surface
    print("Step 3: The shape of the water surface is the geometric cross-section created by this central plane.")
    print("        This plane is parallel to the triangular base face.")

    # Step 4: State the final shape
    shape_name = "regular decagon"
    number_of_sides = 10
    print(f"\nResult:")
    print(f"The cross-section of a regular icosahedron through its center and parallel to a face is a {shape_name}.")
    print(f"This shape is a polygon with {number_of_sides} equal sides and angles.")

    final_answer = f"The shape of the water surface will be a {shape_name}."
    print("\nFinal Answer: " + final_answer)


solve_icosahedron_puzzle()