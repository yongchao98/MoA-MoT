def solve_polytope_f_vector():
    """
    This function identifies and describes the f-vector of the non-simplicial
    4-polytope with 6 vertices and the maximal number of 2-faces.

    According to the classification of 4-polytopes, there is only one
    non-simplicial 4-polytope with 6 vertices: the pyramid over a square pyramid.
    Therefore, it trivially has the maximal number of 2-faces for its class.
    """
    # The f-vector components for the pyramid over a square pyramid.
    # f0: vertices
    # f1: edges
    # f2: 2-faces (polygons)
    # f3: 3-faces (facets/cells)
    f0 = 6
    f1 = 13
    f2 = 13
    f3 = 6

    print("The 4-polytope is the pyramid over a square pyramid.")
    print(f"Its f-vector (f0, f1, f2, f3) is: ({f0}, {f1}, {f2}, {f3}).\n")

    print("To verify, we check the Euler-Poincaré formula for 4-polytopes (f0 - f1 + f2 - f3 = 0):")
    result = f0 - f1 + f2 - f3
    
    # Print the equation with the actual numbers
    print(f"{f0} - {f1} + {f2} - {f3} = {result}")

    if result == 0:
        print("The formula holds true.")
    else:
        print("The formula does not hold, there might be an error.")

solve_polytope_f_vector()