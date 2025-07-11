import sys

def solve_f_vector():
    """
    Calculates the f-vector for the 4-polytope with 6 vertices and the maximal number of faces.
    """
    # The user asks for the f-vector of a "non-simplicial 4-polytope with 6 vertices".
    # However, a key theorem in polytope theory states that any d-polytope with d+2 vertices
    # is necessarily simplicial. For d=4, any 4-polytope with 4+2=6 vertices must be simplicial.
    # This means a non-simplicial 4-polytope with 6 vertices does not exist.
    #
    # The condition "maximal number of 2-faces" points to the Upper Bound Theorem, which is
    # achieved by simplicial neighborly polytopes. Therefore, we interpret the question as asking
    # for the f-vector of the cyclic 4-polytope with 6 vertices, C(6, 4), which is the
    # simplicial polytope that maximizes all face counts for n=6.

    n = 6  # Number of vertices
    d = 4  # Dimension

    # f_0: Number of vertices (0-faces)
    f0 = n

    # f_1: Number of edges (1-faces). For a cyclic/neighborly 4-polytope,
    # every pair of vertices is connected by an edge.
    f1 = n * (n - 1) // 2

    # f_2: Number of 2-faces. The formula for a cyclic 4-polytope is n*(n-3).
    f2 = n * (n - 3)

    # f_3: Number of 3-faces (cells). The formula for a cyclic 4-polytope is n*(n-3)/2.
    f3 = n * (n - 3) // 2

    f_vector = (f0, f1, f2, f3)

    # Verify with Euler's Polyhedral Formula for 4D: f0 - f1 + f2 - f3 = 0
    euler_char = f0 - f1 + f2 - f3

    print(f"A non-simplicial 4-polytope with 6 vertices does not exist.")
    print(f"Assuming the question meant the simplicial polytope with the maximal number of faces (the cyclic polytope C(6, 4)):")
    print("-" * 30)
    print(f"The f-vector (f_0, f_1, f_2, f_3) is: {f_vector}")
    print(f"f_0 (vertices) = {f0}")
    print(f"f_1 (edges) = {f1}")
    print(f"f_2 (2-faces) = {f2}")
    print(f"f_3 (cells) = {f3}")
    print("-" * 30)
    print("Verification using Euler's formula (f_0 - f_1 + f_2 - f_3 = 0):")
    print(f"{f0} - {f1} + {f2} - {f3} = {euler_char}")

    # Final answer in the required format
    # The content is the f-vector tuple.
    sys.stdout.write(f"\n<<<{f_vector}>>>\n")

solve_f_vector()