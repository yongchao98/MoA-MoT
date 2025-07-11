def solve_triangulation():
    """
    Calculates the number of triangles on the outer hull of a great icosahedron.
    
    The outer hull of a great icosahedron is a regular icosahedron. This function
    calculates the number of faces of a regular icosahedron using Euler's formula:
    V - E + F = 2
    """
    
    # A regular icosahedron has 12 vertices.
    V = 12
    
    # At each vertex, 5 edges meet. Each edge connects two vertices.
    # So, the number of edges E = (V * 5) / 2.
    E = (V * 5) / 2
    
    # According to Euler's formula for polyhedra, V - E + F = 2.
    # We can find the number of faces (F) with F = 2 - V + E.
    # Since the faces of an icosahedron are triangles, F is our answer.
    F = 2 - V + E
    
    print("The outer hull of the great icosahedron is a regular icosahedron.")
    print("To find the number of its triangular faces, we use Euler's formula: V - E + F = 2")
    print(f"Number of Vertices (V) = {V}")
    print(f"Number of Edges (E) = ({V} * 5) / 2 = {int(E)}")
    print("\nCalculating the number of faces (F):")
    print("F = 2 - V + E")
    print(f"F = 2 - {V} + {int(E)}")
    final_result = 2 - V + int(E)
    print(f"F = {final_result}")
    
    print(f"\nThe minimal number of triangles needed for triangulating the outer-hull is {final_result}.")

solve_triangulation()
<<<20>>>