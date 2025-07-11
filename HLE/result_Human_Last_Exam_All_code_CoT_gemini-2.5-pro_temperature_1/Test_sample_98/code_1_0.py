def solve_icosahedron_surface_shape():
    """
    This function explains and determines the shape of the water surface
    in a half-filled icosahedron tank standing on one of its faces.
    """

    # Define the reasoning steps
    reasoning = [
        "1. An icosahedron is a regular polyhedron with 20 identical equilateral triangular faces. When it stands on one face, that face is the horizontal base.",
        "2. When half-filled, the water occupies half the total volume. Due to the icosahedron's symmetry, the water surface corresponds to the horizontal plane that cuts the shape at its exact mid-height.",
        "3. This mid-plane is parallel to the top and bottom triangular faces.",
        "4. By analyzing the structure of the icosahedron, this mid-plane is found to intersect six of the polyhedron's side faces.",
        "5. The intersection creates a six-sided polygon (a hexagon).",
        "6. The inherent symmetries of the icosahedron (specifically, 3-fold rotational symmetry around the vertical axis) require this hexagon to be regular, meaning all its sides and angles are equal."
    ]

    # The final answer derived from the reasoning
    final_shape = "A regular hexagon"

    # Print the explanation and the final answer
    print("Reasoning for the shape of the water surface:")
    for step in reasoning:
        print(f"- {step}")
    
    print("\nConclusion:")
    print("The shape of the water surface will be:")
    print(final_shape)

# Execute the function to get the answer
solve_icosahedron_surface_shape()