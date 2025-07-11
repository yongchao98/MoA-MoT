import math

def count_nonzero_terms():
    """
    This function calculates the number of non-zero terms in the asymptotic expansion
    of f(x) up to and including the term in x^-100.
    """
    limit = 100
    
    # We count the number of zero terms first.
    # A term a_n is zero if n belongs to one of two mutually exclusive cases.

    # Case 1: n is an even power of 2 (n = 2^k, k is even).
    # n = 2^(2*j) for j = 0, 1, 2, ...
    zeros_case1_list = []
    j = 0
    while True:
        n = 2**(2 * j)
        if n > limit:
            break
        zeros_case1_list.append(n)
        j += 1
    count_case1 = len(zeros_case1_list)
    
    # Case 2: n = m * 2^k where m > 1 is odd and k is odd.
    zeros_case2_list = []
    k = 1  # Start with the first odd exponent
    while True:
        power_of_2 = 2**k
        # Since m must be at least 3, if 3 * power_of_2 > limit, we can stop.
        if 3 * power_of_2 > limit:
            break
        
        m = 3  # Start with the first odd integer > 1
        while True:
            n = m * power_of_2
            if n > limit:
                break
            zeros_case2_list.append(n)
            m += 2
        
        k += 2  # Move to the next odd exponent
    count_case2 = len(zeros_case2_list)
    
    # The total number of zero terms is the sum of counts from both cases.
    total_zeros = count_case1 + count_case2
    
    # The number of non-zero terms is the total number of terms minus the zero terms.
    non_zero_terms_count = limit - total_zeros
    
    # Output the details of the calculation
    print(f"The number of terms considered is up to n={limit}.")
    print(f"Number of zero terms from Case 1 (n = 2^k, k even): {count_case1}")
    print(f"Number of zero terms from Case 2 (n = m * 2^k, m>1 odd, k odd): {count_case2}")
    print(f"Total number of zero terms is the sum of the two cases:")
    print(f"{count_case1} + {count_case2} = {total_zeros}")
    print(f"The number of non-zero terms is the total number of terms minus the number of zero terms:")
    print(f"{limit} - {total_zeros} = {non_zero_terms_count}")

count_nonzero_terms()