def calculate_euler_characteristic_for_triangulated_surface(V, E, F):
    """
    Calculates the Euler characteristic for a triangulated surface.
    The Euler characteristic is given by the formula chi = V - E + F,
    where V is the number of vertices, E is the number of edges, and
    F is the number of faces.

    Args:
        V (int): The number of vertices.
        E (int): The number of edges.
        F (int): The number of faces.
    """
    chi = V - E + F
    print("The surface is a triangulation of a 2-sphere (a tetrahedron).")
    print(f"Number of Vertices (V): {V}")
    print(f"Number of Edges (E): {E}")
    print(f"Number of Faces (F): {F}")
    print("\nThe Euler characteristic equation is: chi = V - E + F")
    print(f"chi = {V} - {E} + {F} = {chi}")
    print("\nSince the Euler characteristic is non-zero, any map f: S^2 -> S^2")
    print("homotopic to the identity must have a fixed point. This obstructs")
    print("the existence of a section for the fibration conf_2(S^2) -> S^2.")


# A tetrahedron is a simple triangulation of the 2-sphere S^2.
vertices = 4
edges = 6
faces = 4

calculate_euler_characteristic_for_triangulated_surface(vertices, edges, faces)