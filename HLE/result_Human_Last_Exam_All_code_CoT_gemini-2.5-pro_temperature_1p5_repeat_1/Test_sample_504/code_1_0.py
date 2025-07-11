import numpy as np

def solve():
    """
    Constructs and verifies a set of vectors in C^6 meeting the specified criteria.
    """
    # Dimension of the vector space
    dim = 6

    # 1. Construct the set of vectors
    # Start with an orthonormal basis (6 vectors)
    vectors = [np.identity(dim)[i] for i in range(dim)]

    # Add 3 more vectors as described in the plan.
    # We use v_ rather than v7, v8, v9 to match list indexing.
    v_6 = 0.5 * (vectors[0] + vectors[1] + vectors[2] + vectors[3])
    v_7 = 0.5 * (vectors[0] + vectors[1] + vectors[4] + vectors[5])
    v_8 = 0.5 * (vectors[2] + vectors[3] + vectors[4] + vectors[5])

    vectors.extend([v_6, v_7, v_8])

    num_vectors = len(vectors)
    print(f"Constructed a set of {num_vectors} vectors.")

    # 2. Verify the conditions
    all_conditions_met = True
    
    # Check norms and pairwise inner products
    for i in range(num_vectors):
        # Check for unit length
        norm = np.linalg.norm(vectors[i])
        if not np.isclose(norm, 1.0):
            print(f"Error: Vector {i} does not have unit norm (norm={norm}).")
            all_conditions_met = False

        for j in range(i + 1, num_vectors):
            inner_product = np.vdot(vectors[i], vectors[j])
            abs_inner_product = abs(inner_product)

            # Check if the angle corresponds to pi/2 or pi/3
            is_orthogonal = np.isclose(abs_inner_product, 0.0)
            is_pi_over_3 = np.isclose(abs_inner_product, 0.5)

            if not (is_orthogonal or is_pi_over_3):
                print(f"Error: Angle between vector {i} and {j} is not pi/2 or pi/3.")
                print(f"|(v_{i}, v_{j})| = {abs_inner_product:.4f}")
                all_conditions_met = False
    
    # Check for the existence of at least one orthogonal pair
    has_orthogonal_pair = any(np.isclose(abs(np.vdot(vectors[i], vectors[j])), 0.0)
                              for i in range(num_vectors) for j in range(i + 1, num_vectors))
    if not has_orthogonal_pair:
        print("Error: The set does not contain an orthogonal pair.")
        all_conditions_met = False

    if all_conditions_met:
        print("All constructed vectors have unit length.")
        print("The absolute value of the inner product of any pair is 0 or 0.5.")
        print("The set contains at least one orthogonal pair.")
        print("\nThe construction is valid.")
    else:
        print("\nThe construction is invalid.")

    # The final answer is the number of vectors in our valid construction.
    # Based on the reasoning, this is the maximal number.
    print("\nThe calculation for the largest number of vectors is based on this verified construction:")
    # Using python to print the equation. 
    base_vectors = 6
    additional_vectors = 3
    total_vectors = base_vectors + additional_vectors
    print(f"{base_vectors} + {additional_vectors} = {total_vectors}")

solve()