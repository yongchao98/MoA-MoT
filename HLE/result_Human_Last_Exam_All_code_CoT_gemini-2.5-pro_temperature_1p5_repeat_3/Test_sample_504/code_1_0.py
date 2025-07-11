import numpy as np

def solve():
    """
    This function constructs and verifies a set of vectors in C^6 to find the largest
    number of such vectors where the angle between any two is pi/2 or pi/3,
    and at least one orthogonal pair exists.
    """
    # Step 1: Define the 9 vectors of the Hesse SIC-POVM in C^3.
    w = np.exp(2j * np.pi / 3)
    # The 9 vectors are the columns of this 3x9 matrix, scaled by 1/sqrt(2).
    sic_matrix_T = np.array([
        [1, w, w**2, 0, 0, 0, 1, w, w**2],
        [0, 0, 0, 1, w, w**2, 1, w**2, w],
        [1, w**2, w, 1, w**2, w, 0, 0, 0]
    ], dtype=np.complex128)
    
    U = (1 / np.sqrt(2)) * sic_matrix_T.T
    W = U  # We can use the same set for the second C^3 subspace.

    # Step 2: Define the three special vectors in C^3.
    # These are derived from the standard basis vectors.
    A = (1 / np.sqrt(2)) * np.eye(3, dtype=np.complex128)
    B = A # Use the same set for the second C^3 subspace.
    
    num_u = U.shape[0]
    num_a = A.shape[0]

    # Step 3: Construct the 27 vectors in C^6.
    vectors = []
    
    # S1: 9 vectors of the form (u, 0)
    for i in range(num_u):
        v = np.zeros(6, dtype=np.complex128)
        v[:3] = U[i]
        vectors.append(v)
    
    # S2: 9 vectors of the form (0, w)
    for i in range(num_u):
        v = np.zeros(6, dtype=np.complex128)
        v[3:] = W[i]
        vectors.append(v)

    # S3: 9 vectors of the form (a, b)
    for i in range(num_a):
        for j in range(num_a):
            v = np.zeros(6, dtype=np.complex128)
            v[:3] = A[i]
            v[3:] = B[j]
            vectors.append(v)

    # Step 4: Verify the properties of the constructed set.
    n = len(vectors)
    has_orthogonal_pair = False
    all_conditions_met = True
    
    for i in range(n):
        for j in range(i, n):
            # Inner product is conj(v1) dot v2
            inner_product = np.dot(np.conj(vectors[i]), vectors[j])
            abs_inner_product = np.abs(inner_product)
            
            if i == j: # Check normalization
                if not np.isclose(abs_inner_product, 1.0):
                    print(f"Vector {i} is not normalized: |v|^2 = {abs_inner_product**2}")
                    all_conditions_met = False
            else: # Check angle conditions
                is_pi_2 = np.isclose(abs_inner_product, 0.0)
                is_pi_3 = np.isclose(abs_inner_product, 0.5)
                
                if not (is_pi_2 or is_pi_3):
                    print(f"Invalid angle between vector {i} and {j}: |(v_i, v_j)| = {abs_inner_product}")
                    all_conditions_met = False
                
                if is_pi_2:
                    has_orthogonal_pair = True

    if not has_orthogonal_pair:
        print("The set does not contain an orthogonal pair.")
        all_conditions_met = False

    if all_conditions_met:
        print("The construction is successful and all conditions are met.")
        print("The total number of vectors is the sum of vectors from three sets:")
        num_s1 = num_u
        num_s2 = num_u
        num_s3 = num_a * num_a
        print(f"{num_s1} + {num_s2} + {num_s3} = {n}")
        print("\nThe largest number of such vectors is:")
        print(n)
    else:
        print("\nThe construction failed the verification.")

solve()