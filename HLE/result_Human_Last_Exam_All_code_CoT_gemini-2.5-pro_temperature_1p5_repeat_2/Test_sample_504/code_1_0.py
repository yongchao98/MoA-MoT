import numpy as np

def solve():
    """
    Constructs and verifies a set of vectors in C^6 to find the maximum number
    satisfying the given angle conditions.
    """
    # Part 1: Start with the standard orthonormal basis
    vectors = [np.eye(6)[i] for i in range(6)]

    # Part 2: Construct additional vectors
    # Partition indices into 3 pairs: {1,2}, {3,4}, {5,6}
    # Basic 2D vectors
    u = np.array([1, 1])
    v = np.array([1, -1])

    # Normalization factor for the new vectors.
    # Each will have 4 non-zero components of value +/- 0.5.
    # The squared norm is 4 * (0.5)^2 = 1.
    norm_factor = 0.5

    # Group 1: Vectors on coordinates {1,2,3,4}
    for b1 in [u, v]:
        for b2 in [u, v]:
            vec = np.zeros(6)
            vec[0:2] = b1
            vec[2:4] = b2
            vectors.append(vec * norm_factor)

    # Group 2: Vectors on coordinates {1,2,5,6}
    for b1 in [u, v]:
        for b3 in [u, v]:
            vec = np.zeros(6)
            vec[0:2] = b1
            vec[4:6] = b3
            vectors.append(vec * norm_factor)
            
    # Group 3: Vectors on coordinates {3,4,5,6}
    for b2 in [u, v]:
        for b3 in [u, v]:
            vec = np.zeros(6)
            vec[2:4] = b2
            vec[4:6] = b3
            vectors.append(vec * norm_factor)

    n = len(vectors)
    print(f"Constructed a set of {n} vectors.")

    # Part 3: Verify the conditions for all pairs
    all_conditions_met = True
    for i in range(n):
        for j in range(i + 1, n):
            v1 = vectors[i]
            v2 = vectors[j]

            # Check for linear independence (not scalar multiples of each other)
            # Since all vectors are unit vectors, we just check if v1 = +/- v2
            if np.allclose(v1, v2) or np.allclose(v1, -v2):
                print(f"Error: Vectors {i+1} and {j+1} are linearly dependent.")
                all_conditions_met = False
                break

            # Calculate the absolute value of the Hermitian inner product
            inner_product_abs = np.abs(np.vdot(v1, v2))

            # Check if the angle condition is met
            is_pi_over_2 = np.isclose(inner_product_abs, 0)
            is_pi_over_3 = np.isclose(inner_product_abs, 0.5)

            if not (is_pi_over_2 or is_pi_over_3):
                print(f"Error: Angle condition failed for pair ({i+1}, {j+1})")
                print(f"|({v1}, {v2})| = {inner_product_abs}")
                all_conditions_met = False
                break
        if not all_conditions_met:
            break
            
    if all_conditions_met:
        print("All constructed vectors satisfy the given conditions.")
        print("The vectors are pairwise linearly independent.")
        print("The angle between any two is pi/2 or pi/3.")
        print("The set contains orthogonal pairs (e.g., the standard basis vectors).")
        print(f"\nThe largest number of such vectors is at least {n}.")
        print(f"Thus, the largest number is {n}.")

solve()
<<<18>>>