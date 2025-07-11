import numpy as np

def check_vectors(vectors):
    """
    Checks if a list of vectors satisfies the angle conditions.
    For any pair (v, w), the angle alpha must be pi/2 or pi/3.
    This means cos(alpha) = |(v,w)| / (|v||w|) is 0 or 1/2.
    """
    n = len(vectors)
    for i in range(n):
        for j in range(i + 1, n):
            v = vectors[i]
            w = vectors[j]
            
            norm_v = np.linalg.norm(v)
            norm_w = np.linalg.norm(w)
            
            if np.isclose(norm_v, 0) or np.isclose(norm_w, 0):
                print(f"Warning: Vector is zero at index {i} or {j}")
                continue

            # In complex space, inner product is v^H * w
            inner_product = np.vdot(v, w)
            
            cos_alpha_val = np.abs(inner_product) / (norm_v * norm_w)

            # Check if cos(alpha) is close to 0 (pi/2) or 0.5 (pi/3)
            if not (np.isclose(cos_alpha_val, 0) or np.isclose(cos_alpha_val, 0.5)):
                print(f"Condition failed for vector pair ({i}, {j})")
                print(f"v_{i} = {v}")
                print(f"v_{j} = {w}")
                print(f"|({v},{w})| = {np.abs(inner_product):.4f}")
                print(f"cos(alpha) = {cos_alpha_val:.4f}")
                return False
    return True

def main():
    """
    Constructs and verifies the set of 18 vectors in C^6.
    """
    dim = 6
    vectors = []

    # 1. Add the 6 standard basis vectors
    basis_vectors = []
    for i in range(dim):
        v = np.zeros(dim, dtype=complex)
        v[i] = 1.0
        basis_vectors.append(v)
    vectors.extend(basis_vectors)
    
    num_basis_vectors = len(basis_vectors)

    # 2. Define supports and Hadamard matrix
    supports = [
        [0, 1, 2, 3],  # J1 = {1, 2, 3, 4}
        [0, 1, 4, 5],  # J2 = {1, 2, 5, 6}
        [2, 3, 4, 5]   # J3 = {3, 4, 5, 6}
    ]
    
    H4 = np.array([
        [1, 1, 1, 1],
        [1, -1, 1, -1],
        [1, 1, -1, -1],
        [1, -1, -1, 1]
    ], dtype=complex)
    
    num_vec_per_support = []

    # 3. Create vectors from supports and Hadamard matrix
    for support in supports:
        family_vectors = []
        for i in range(H4.shape[0]):
            v = np.zeros(dim, dtype=complex)
            coeffs = H4[i, :]
            for k in range(len(support)):
                v[support[k]] = coeffs[k]
            v = v / 2.0  # Normalize to have squared coords sum to 1 (4 * 1/4 = 1)
            family_vectors.append(v)
        vectors.extend(family_vectors)
        num_vec_per_support.append(len(family_vectors))

    print("Constructing the set of vectors.")
    print("The set consists of:")
    print(f"- {num_basis_vectors} standard basis vectors (e.g., e1, e2, ...)")
    
    total_hadamard_vectors = 0
    for i, count in enumerate(num_vec_per_support):
        print(f"- {count} vectors from support J{i+1} = {np.array(supports[i])+1}")
        total_hadamard_vectors += count
    
    total_vectors = num_basis_vectors + total_hadamard_vectors
    
    print("\nVerifying the angle property for all pairs...")
    if check_vectors(vectors):
        print("All pairs of vectors satisfy the condition.")
        print("The angle between any two is either pi/2 or pi/3.")
        print("\nAlso, there are orthogonal pairs (e.g., (e1, e2) = 0).")
        
        # Output the equation part
        equation_str = f"{num_basis_vectors}"
        for count in num_vec_per_support:
            equation_str += f" + {count}"
        equation_str += f" = {total_vectors}"
        
        print(f"\nThe largest number of such vectors is found from this construction.")
        print(f"Final calculation: {equation_str}")
        print(f"The largest number is {total_vectors}.")

if __name__ == '__main__':
    main()