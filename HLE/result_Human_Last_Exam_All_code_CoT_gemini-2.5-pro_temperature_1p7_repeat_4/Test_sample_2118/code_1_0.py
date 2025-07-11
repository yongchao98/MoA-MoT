def v2(n):
    """
    Calculates the exponent of 2 in the prime factorization of n (2-adic valuation).
    """
    if n == 0:
        return -1 # Or float('inf'), an arbitrary choice for this case
    if n % 2 != 0:
        return 0
    count = 0
    # In CPython, this loop is efficient for powers of 2.
    while n > 0 and n % 2 == 0:
        n //= 2
        count += 1
    return count

def count_nonzero_terms():
    """
    Counts the number of integers n in [2, 100] for which v2(n) is even.
    This corresponds to the number of non-zero coefficients c_n.
    """
    limit = 100
    counts = {} # To store counts for each even exponent value

    for n in range(2, limit + 1):
        exponent = v2(n)
        if exponent % 2 == 0:
            if exponent not in counts:
                counts[exponent] = 0
            counts[exponent] += 1
    
    total_count = 0
    sum_parts = []
    
    # Sort by exponent to present the result in a structured way
    sorted_exponents = sorted(counts.keys())
    
    print("The number of nonzero terms is determined by counting integers n from 2 to 100")
    print("where the exponent of 2 in their prime factorization, v2(n), is even.")
    print("-" * 30)

    for exponent in sorted_exponents:
        count = counts[exponent]
        if exponent == 0:
            print(f"Number of terms where n is odd (v2(n) = 0): {count}")
        else:
            print(f"Number of terms where v2(n) = {exponent}: {count}")
        total_count += count
        sum_parts.append(str(count))
        
    final_equation = " + ".join(sum_parts)
    print("-" * 30)
    print(f"The total number of nonzero terms is the sum: {final_equation} = {total_count}")

if __name__ == '__main__':
    count_nonzero_terms()
