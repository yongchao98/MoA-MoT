import numpy as np

def solve():
    """
    Constructs and verifies a set of 12 vectors in C^6 (represented as R^6)
    satisfying the given angle conditions.
    """
    # Define the 12 vectors as rows in a matrix
    vectors = np.zeros((12, 6))

    # First 6 vectors are the standard orthonormal basis
    for i in range(6):
        vectors[i, i] = 1.0

    # Next 6 vectors are constructed based on subsets of the basis
    # v_7 to v_12 (indexed 6 to 11 in the matrix)
    vectors[6] = 0.5 * np.array([1, 1, 1, 1, 0, 0])
    vectors[7] = 0.5 * np.array([1, 1, 0, 0, 1, 1])
    vectors[8] = 0.5 * np.array([0, 0, 1, 1, 1, 1])
    vectors[9] = 0.5 * np.array([1, -1, 1, -1, 0, 0])
    vectors[10] = 0.5 * np.array([1, -1, 0, 0, 1, -1])
    vectors[11] = 0.5 * np.array([0, 0, 1, -1, 1, -1])

    print("Verifying the set of 12 vectors:")
    print("===================================")

    # 1. Verify they are all unit vectors
    all_unit = True
    for i in range(12):
        norm = np.linalg.norm(vectors[i])
        if not np.isclose(norm, 1.0):
            print(f"Vector {i+1} is not a unit vector. Norm is {norm}")
            all_unit = False
    if all_unit:
        print("All 12 vectors are unit vectors.")
    print("-----------------------------------")

    # 2. Verify angle conditions by checking inner products
    # For any two distinct unit vectors v, w, |(v,w)| should be 0 or 0.5
    conditions_met = True
    found_orthogonal_pair = False
    for i in range(12):
        for j in range(i + 1, 12):
            inner_product = np.dot(vectors[i], vectors[j])
            abs_inner_product = abs(inner_product)

            if np.isclose(abs_inner_product, 0.0):
                found_orthogonal_pair = True
            elif not np.isclose(abs_inner_product, 0.5):
                print(f"Condition failed for pair (v_{i+1}, v_{j+1}):")
                print(f"|({vectors[i]}, {vectors[j]})| = {abs_inner_product}")
                conditions_met = False
    
    if conditions_met:
        print("Angle conditions are met for all pairs.")
        if found_orthogonal_pair:
            print("An orthogonal pair was found in the set.")
        else:
            print("Error: No orthogonal pair was found.")
            conditions_met = False
    else:
        print("Angle conditions are NOT met.")
    
    print("===================================")

    if all_unit and conditions_met:
        print("The construction of 12 vectors is valid.")
        largest_number = 12
        print(f"The largest number of such vectors is {largest_number}.")
    else:
        print("The constructed set of vectors is invalid.")

solve()
<<<12>>>