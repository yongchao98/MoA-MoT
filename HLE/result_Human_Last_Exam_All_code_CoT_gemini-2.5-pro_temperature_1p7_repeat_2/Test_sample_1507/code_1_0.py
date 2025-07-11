def solve_great_icosahedron_triangulation():
    """
    Calculates the number of triangles on the visible outer surface of the great icosahedron.
    """
    
    # The visible surface of the great icosahedron is a stellation of the regular icosahedron.
    # It is formed by placing a triangular pyramid on each face of a central icosahedron.
    
    # 1. A regular icosahedron has 20 triangular faces.
    # This is the number of pyramids that form the visible surface.
    icosahedron_faces = 20
    
    # 2. Each spike is a triangular pyramid. One of its faces serves as the base,
    # and the other three faces are the visible sides.
    visible_faces_per_pyramid = 3
    
    # 3. The total number of visible triangles is the product of the number of 
    # icosahedron faces and the number of visible faces per pyramid.
    # Since the surface is already made of triangles, this is the minimal number
    # for its triangulation.
    total_visible_triangles = icosahedron_faces * visible_faces_per_pyramid
    
    print("To find the minimal number of triangles for triangulating the outer-hull of the great icosahedron:")
    print(f"1. We start with a central icosahedron, which has {icosahedron_faces} faces.")
    print(f"2. A pyramid is built on each face, and each pyramid has {visible_faces_per_pyramid} visible triangular faces.")
    print("3. The total number of visible triangles is the product of these two numbers.")
    print("\nFinal Calculation:")
    print(f"{icosahedron_faces} * {visible_faces_per_pyramid} = {total_visible_triangles}")

solve_great_icosahedron_triangulation()
<<<60>>>