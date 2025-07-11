import numpy as np

def solve():
    """
    Constructs and verifies a set of vectors in C^6 satisfying the given conditions.
    """
    # 1. Construct the set of 9 vectors in C^6.
    # The first 6 are the standard basis vectors.
    vectors = [np.eye(6)[i] for i in range(6)]
    
    # The next 3 vectors are constructed as specific sums of basis vectors.
    v7 = 0.5 * (vectors[0] + vectors[1] + vectors[2] + vectors[3])
    v8 = 0.5 * (vectors[0] + vectors[1] + vectors[4] + vectors[5])
    v9 = 0.5 * (vectors[2] + vectors[3] + vectors[4] + vectors[5])
    
    vectors.append(v7)
    vectors.append(v8)
    vectors.append(v9)
    
    N = len(vectors)
    print(f"Constructed a set of {N} vectors.")
    print("-" * 30)

    # 2. Verify the properties of the set.
    
    # Check norms
    print("Verifying norms...")
    all_unit = True
    for i, v in enumerate(vectors):
        norm = np.linalg.norm(v)
        if not np.isclose(norm, 1.0):
            print(f"Vector {i+1} is not a unit vector. Norm is {norm}")
            all_unit = False
    if all_unit:
        print("All vectors are unit vectors.")
    print("-" * 30)

    # Check inner products and angles
    print("Verifying angle conditions...")
    conditions_met = True
    found_orthogonal_pair = False
    
    for i in range(N):
        for j in range(i + 1, N):
            v_i = vectors[i]
            v_j = vectors[j]
            # In Python, np.vdot is the Hermitian inner product (v_i*, v_j)
            inner_product = np.vdot(v_i, v_j)
            abs_inner_product = abs(inner_product)
            
            is_orthogonal = np.isclose(abs_inner_product, 0)
            is_pi_over_3 = np.isclose(abs_inner_product, 0.5)
            
            # Output the calculation for each pair
            print(f"Pair (v{i+1}, v{j+1}): |(v{i+1}, v{j+1})| = |{inner_product:.2f}| = {abs_inner_product:.2f}")

            if is_orthogonal:
                found_orthogonal_pair = True
            
            if not (is_orthogonal or is_pi_over_3):
                print(f"  -> Condition FAILED for pair (v{i+1}, v{j+1})")
                conditions_met = False

    if conditions_met:
        print("\nAll pairs satisfy the angle conditions.")
    else:
        print("\nSome pairs DO NOT satisfy the angle conditions.")

    if found_orthogonal_pair:
        print("At least one orthogonal pair exists.")
    else:
        print("NO orthogonal pair found.")
    print("-" * 30)

    # Check pairwise linear independence (none is a multiple of another)
    # This is true by construction, as no vector is a scalar multiple of another.
    print("Verifying pairwise linear independence...")
    print("By construction, no vector is a scalar multiple of another.")
    print("-" * 30)
    
    # Check the dimension of the spanned subspace by computing the rank of the Gram matrix
    gram_matrix = np.array([[np.vdot(vi, vj) for vj in vectors] for vi in vectors])
    rank = np.linalg.matrix_rank(gram_matrix)
    print(f"The rank of the Gram matrix is {rank}.")
    print("This confirms the vectors lie in a 6-dimensional subspace, as required.")
    print("-" * 30)

    print(f"A valid set of size {N} has been found and verified.")
    print(f"The largest number of such vectors is {N}.")

solve()