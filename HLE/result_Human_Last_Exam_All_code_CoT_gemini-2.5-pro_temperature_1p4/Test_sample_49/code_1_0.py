import math

def solve_f_vector():
    """
    Calculates the f-vector of the non-simplicial 4-polytope with 6 vertices
    and the maximal number of 2-faces.
    """
    # Step 1: We are given the number of vertices.
    f0 = 6

    # Step 2: The maximal number of edges for a polytope with 6 vertices is C(6,2).
    # To maximize the number of 2-faces, the polytope should have the maximal number of edges.
    f1 = int(math.comb(6, 2))

    # Step 3: For any 4-polytope with 6 vertices, the absolute maximum number of 2-faces is 18,
    # achieved by the simplicial neighborly polytope with f-vector (6, 15, 18, 9).
    # Since we are looking for a non-simplicial polytope, its number of 2-faces must be
    # strictly less than 18. The maximal integer value is therefore 17.
    f2 = 18 - 1

    # Step 4: Use Euler's formula for 4-polytopes to find f3.
    # The formula is: f0 - f1 + f2 - f3 = 0
    # So, f3 = f0 - f1 + f2
    f3 = f0 - f1 + f2

    # Step 5: Print the results and the final equation.
    print("The f-vector of the non-simplicial 4-polytope with 6 vertices and the maximal number of 2-faces is (f0, f1, f2, f3).")
    print(f"\nGiven values and deductions:")
    print(f"f0 (vertices) = {f0}")
    print(f"f1 (edges) = {f1} (maximal possible for 6 vertices)")
    print(f"f2 (2-faces) = {f2} (maximal for a non-simplicial case)")

    print("\nUsing Euler's formula for 4-polytopes to find f3:")
    print(f"f0 - f1 + f2 - f3 = 0")
    print(f"{f0} - {f1} + {f2} - f3 = 0")
    sum_val = f0 - f1 + f2
    print(f"{sum_val} - f3 = 0")
    print(f"f3 = {f3}")

    print(f"\nThe final f-vector is: ({f0}, {f1}, {f2}, {f3})")

solve_f_vector()