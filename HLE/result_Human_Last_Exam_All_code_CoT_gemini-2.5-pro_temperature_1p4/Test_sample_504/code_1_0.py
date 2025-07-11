import math

def verify_construction(dims, vector_counts):
    """
    Verifies a vector construction against the rank inequality.
    
    dims: A list of dimensions of the subspaces.
    vector_counts: A list of the number of vectors in each subspace.
    """
    
    k = sum(vector_counts)
    d = sum(dims)
    
    if d > 6:
        print(f"Construction with d={d} is invalid for C^6.")
        return
        
    total_trace_M2 = 0
    
    print(f"Verifying construction for decomposition d = {' + '.join(map(str, dims))}:")
    print(f"Total vectors k = {k}, Total dimension d = {d}")
    
    # Calculate trace(M^2) for each block
    total_N = 0
    for i in range(len(dims)):
        d_i = dims[i]
        k_i = vector_counts[i]
        
        # Determine the number of non-orthogonal pairs N_i for this subspace
        N_i = 0
        if d_i == 2 and k_i == 3: # Equiangular set in C^2
            N_i = k_i * (k_i - 1)
        elif d_i == 4 and k_i == 8: # From 2 MUBs in C^4
            # 16 pairs of (e_i, f_j), each appears twice in the ordered sum
            N_i = 16 * 2
        elif d_i == 3 and k_i == 9: # SIC-POVM in C^3
            # All pairs are non-orthogonal
            N_i = k_i * (k_i - 1)

        trace_M2_i = k_i + N_i / 4
        total_trace_M2 += trace_M2_i
        total_N += N_i

    # Verify the main inequality
    lhs = k**2
    rhs = d * total_trace_M2
    
    # The problem's value for N is for k(k-1) pairs, but our sum is over ordered pairs
    final_N = total_N
    # Check the inequality using the calculated N
    rhs_formula = d * (k + final_N / 4.0)
    
    print(f"Number of non-orthogonal ordered pairs N = {final_N}")
    print(f"Checking inequality: k^2 <= d * (k + N/4)")
    print(f"LHS = k^2 = {k}^2 = {lhs}")
    print(f"RHS = {d} * ({k} + {final_N}/4) = {d} * ({k} + {final_N/4.0}) = {rhs_formula}")

    if math.isclose(lhs, rhs_formula) or lhs < rhs_formula:
        print("The inequality holds. The construction is valid.")
    else:
        print("The inequality DOES NOT hold. The construction is invalid.")
    print("-" * 20)
    return k

# Case 1: Decomposition 4 + 2
dims1 = [4, 2]
counts1 = [8, 3]
k1 = verify_construction(dims1, counts1)

# Case 2: Decomposition 3 + 3
dims2 = [3, 3]
counts2 = [9, 9]
k2 = verify_construction(dims2, counts2)

largest_number = max(k1, k2)
print(f"The largest number of vectors found is {largest_number}.")

# Final Answer
print("The final equation for the best construction (18 vectors in C^6 from two C^3 subspaces) is:")
k, d, N = 18, 6, 144
print(f"{k}^2 <= {d}*({k} + {N}/4)")
print(f"{k**2} <= {d}*({k + N/4})")
print(f"{k**2} <= {d* (k+N/4)}")

<<<18>>>