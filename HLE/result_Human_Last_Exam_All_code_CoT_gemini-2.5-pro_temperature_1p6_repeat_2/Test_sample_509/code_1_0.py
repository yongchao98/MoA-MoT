def calculate_euler_characteristic(manifold_name, vertices, edges, faces):
    """
    Calculates the Euler characteristic for a manifold from its triangulation.
    Formula: chi = V - E + F
    V = number of vertices
    E = number of edges
    F = number of faces
    """
    chi = vertices - edges + faces
    print(f"Calculating Euler characteristic for the {manifold_name}:")
    print(f"Formula: chi = V - E + F")
    # Print the equation with the numbers
    print(f"{vertices} - {edges} + {faces} = {chi}")
    
    if chi == 0:
        print(f"Result: The Euler characteristic is 0. This is consistent with the existence of a non-vanishing vector field, which allows for a section in the configuration space fibration.")
    else:
        print(f"Result: The Euler characteristic is {chi} (non-zero). A section may not exist. For example, any vector field on the sphere must have a zero.")
    print("-" * 20)

def main():
    # Example 1: The Sphere (S^2)
    # We can triangulate the sphere as an octahedron.
    # It has 6 vertices (North pole, South pole, 4 on the equator).
    # It has 12 edges (4 from N to equator, 4 from S to equator, 4 around the equator).
    # It has 8 triangular faces.
    sphere_V = 6
    sphere_E = 12
    sphere_F = 8
    calculate_euler_characteristic("Sphere", sphere_V, sphere_E, sphere_F)
    
    # Example 2: The Torus (T^2)
    # We can triangulate the torus from a square with opposite sides identified.
    # A simple triangulation:
    # It has 1 vertex (all 4 corners of the square are identified).
    # It has 3 edges (the two identified sides 'a' and 'b', and the diagonal 'c').
    # It has 2 faces (the two triangles formed by the diagonal).
    torus_V = 1
    torus_E = 3
    torus_F = 2
    calculate_euler_characteristic("Torus", torus_V, torus_E, torus_F)

if __name__ == "__main__":
    main()
