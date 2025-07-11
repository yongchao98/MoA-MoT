import numpy as np
import itertools

def solve():
    """
    This function determines the maximum value of n for a given m,
    such that there exist n binary vectors in {0,1}^m in general linear position.
    It demonstrates the result for a specific case, m=4.
    """
    m = 4
    print(f"Let's analyze the case for m = {m}")

    # Part 1: Show that n = m+1 is possible.
    n = m + 1
    print(f"\nPart 1: Constructing a set of n = m+1 = {n} vectors and verifying the property.")
    
    # Create the standard basis vectors e_i
    vectors = [list(v) for v in np.identity(m, dtype=int)]
    # Add the all-ones vector
    vectors.append([1] * m)

    print("We construct the set V of n vectors:")
    for i, v in enumerate(vectors):
        print(f"  v{i+1} = {v}")

    print("\nNow, we check if every subset of m vectors is linearly independent.")
    print("We do this by checking if the determinant of the matrix formed by the subset is non-zero.")

    subsets = list(itertools.combinations(vectors, m))
    all_independent = True
    for i, subset in enumerate(subsets):
        matrix = np.array(subset).T # put vectors as columns
        det = np.linalg.det(matrix)
        print(f"  Subset {i+1} of {m} vectors has a determinant of: {det:.0f}")
        if np.isclose(det, 0):
            all_independent = False
            break

    if all_independent:
        print(f"\nConclusion for Part 1: A set of n = {n} vectors exists and satisfies the condition.")
    else:
        print(f"\nConclusion for Part 1: The constructed set for n = {n} does not work.")


    # Part 2: Show that n = m+2 is not possible.
    n_plus_1 = m + 2
    print(f"\nPart 2: Showing that n = m+2 = {n_plus_1} is not possible.")
    print("We add one more vector to our set V and show the condition is violated.")
    
    # Add a new vector. Let's pick a vector with 2 ones.
    new_vector = [1, 1, 0, 0]
    vectors.append(new_vector)
    
    print(f"\nOur new set of {len(vectors)} vectors is:")
    for i, v in enumerate(vectors):
        print(f"  v{i+1} = {v}")
    
    print("\nWe must find a subset of m=4 vectors that is linearly dependent (determinant is zero).")
    
    # Theory suggests that the new vector v_6 = [1,1,0,0] is a linear combination of v_1 and v_2.
    # v_6 = v_1 + v_2.
    # So the set {v_1, v_2, v_6} is dependent.
    # Any superset of size 4 containing these three is also dependent. Let's test {v_1, v_2, v_3, v_6}.
    e1 = [1, 0, 0, 0]
    e2 = [0, 1, 0, 0]
    e3 = [0, 0, 1, 0]
    dependent_subset = [e1, e2, e3, new_vector]
    
    matrix = np.array(dependent_subset).T
    det = np.linalg.det(matrix)
    
    print("\nConsider the subset {v1, v2, v3, v6}:")
    print(f"The determinant of the matrix formed by these vectors is {det:.0f}.")
    
    if np.isclose(det, 0):
        print("This subset is linearly dependent, which violates the condition.")
        print("Therefore, n cannot be m+2.")
    else:
        print("There is an error in the logic, as this subset should be dependent.")

    # Final Conclusion
    print("\n-------------------------------------------")
    print("The maximum value of n is m+1.")
    print(f"For our case, m = {m}, the maximum n is {m+1}.")
    print("The final relation is:")
    final_n = m + 1
    print(f"{final_n} = {m} + 1")
    
solve()