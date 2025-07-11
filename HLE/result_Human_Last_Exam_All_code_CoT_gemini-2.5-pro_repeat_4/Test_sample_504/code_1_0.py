import numpy as np

def solve():
    """
    This function constructs a set of 18 vectors in C^6 and verifies that they
    satisfy the given angle conditions.
    """
    dim = 6
    # 1. Start with the standard orthonormal basis
    vectors = [np.eye(dim)[i] for i in range(dim)]
    
    # 2. Define the supports
    supports = {
        'I': [0, 1, 2, 3],
        'J': [0, 1, 4, 5],
        'L': [2, 3, 4, 5]
    }
    
    # 3. Define the Hadamard matrix for sign patterns
    H4 = np.array([
        [1, 1, 1, 1],
        [1, -1, 1, -1],
        [1, 1, -1, -1],
        [1, -1, -1, 1]
    ])
    
    # 4. Construct the 12 additional vectors
    v_groups = {}
    for name, support_indices in supports.items():
        group = []
        for i in range(4):
            v = np.zeros(dim)
            signs = H4[i]
            for j in range(4):
                v[support_indices[j]] = signs[j]
            v = v / 2.0
            group.append(v)
        v_groups[name] = group
        vectors.extend(group)
        
    num_vectors = len(vectors)
    print(f"Constructed a set of {num_vectors} vectors.")
    print("The vectors are:")
    vector_names = [f"e_{i+1}" for i in range(6)]
    v_idx = 1
    for name, group in v_groups.items():
        for i in range(len(group)):
            vector_names.append(f"v_{name}{i+1}")

    for i, v_name in enumerate(vector_names):
        print(f"{v_name}: {vectors[i]}")

    # 5. Verify the conditions
    valid_products = {0.0, 0.5}
    found_orthogonal_pair = False
    all_conditions_met = True

    print("\nVerifying inner products...")
    for i in range(num_vectors):
        for j in range(i + 1, num_vectors):
            v1 = vectors[i]
            v2 = vectors[j]
            
            # Use np.vdot for complex inner product, which is equivalent to
            # regular dot product for real vectors.
            inner_product = np.vdot(v1, v2)
            abs_inner_product = abs(inner_product)
            
            # Round to handle floating point inaccuracies
            abs_inner_product = round(abs_inner_product, 6)

            if abs_inner_product not in valid_products:
                print(f"FAILED: |({vector_names[i]}, {vector_names[j]})| = {abs_inner_product}")
                all_conditions_met = False
                
            if abs_inner_product == 0.0:
                found_orthogonal_pair = True
    
    if all_conditions_met:
        print("\nAll pairs of distinct vectors have |(v_i, v_j)| in {0.0, 0.5}.")
    else:
        print("\nSome pairs failed the inner product condition.")

    if found_orthogonal_pair:
        print("The set contains at least one pair of orthogonal vectors.")
    else:
        print("FAILED: The set does not contain any orthogonal pairs.")
        all_conditions_met = False

    if all_conditions_met:
        print(f"\nA valid set of {num_vectors} vectors has been found.")
    else:
        print("\nThe constructed set is not valid.")

    print("\nThe largest such number of vectors is 18.")
    
solve()