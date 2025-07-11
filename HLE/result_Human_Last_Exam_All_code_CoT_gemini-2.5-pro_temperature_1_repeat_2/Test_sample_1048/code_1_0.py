import collections

def find_best_permutation(a, x):
    """
    Finds a permutation of 'a' to maximize the final value of 'x' after
    a series of modulo operations.

    Args:
        a (list[int]): A list of n numbers.
        x (int): The initial number.
    """
    if not a:
        print(f"The list is empty. Final x is {x}")
        return

    # 1. Pre-processing: handle numbers larger than current x
    # They don't change x's value, so we can apply them first conceptually
    # and work with a smaller effective list.
    initial_x = x
    
    effective_a = sorted([val for val in a if val <= x])
    if not effective_a:
        # All numbers in 'a' are greater than x, so x never changes.
        # The permutation doesn't matter.
        p = sorted(a)
        print(f"Initial x = {x}")
        for val in p:
            print(f"{x} % {val} = {x % val}")
            x %= val
        print(f"Final x = {x}")
        return

    # From this point, we only need to consider numbers <= initial_x
    a = sorted(list(set(effective_a)))
    n = len(a)
    
    # 2. Find minimum element `m` and the rest of the list `A`
    m = a[0]
    A = a[1:]
    
    # 3. DP state exploration (BFS-style)
    # memo[mask] = set of reachable values using subset 'mask' of A
    memo = {0: {initial_x}}
    # parent[(mask, value)] = (prev_mask, prev_value, number_used) for backtracking
    parent = {}
    
    q = collections.deque([(0, initial_x)])
    
    # Keep track of the state that gives the best remainder modulo m
    max_rem_info = {'rem': initial_x % m, 'y': initial_x, 'mask': 0}

    visited_states = {(0, initial_x)}

    while q:
        mask, y = q.popleft()
        
        for i in range(len(A)):
            if not ((mask >> i) & 1):  # If i-th element of A is not used
                new_mask = mask | (1 << i)
                new_y = y % A[i]
                
                if (new_mask, new_y) not in visited_states:
                    if new_mask not in memo:
                        memo[new_mask] = set()
                    memo[new_mask].add(new_y)
                    parent[(new_mask, new_y)] = (mask, y, A[i])
                    
                    q.append((new_mask, new_y))
                    visited_states.add((new_mask, new_y))

                    rem = new_y % m
                    if rem > max_rem_info['rem']:
                        max_rem_info = {'rem': rem, 'y': new_y, 'mask': new_mask}
    
    # 4. Backtrack to find the best permutation prefix
    p_before = []
    curr_mask = max_rem_info['mask']
    curr_y = max_rem_info['y']
    
    while curr_mask > 0:
        prev_mask, prev_y, val = parent[(curr_mask, curr_y)]
        p_before.insert(0, val)
        curr_mask, curr_y = prev_mask, prev_y
        
    # 5. Construct the full optimal permutation
    p_before_set = set(p_before)
    # Get numbers from original list `a` that were not used in p_before and are not `m`
    all_other_nums = [val for val in a if val != m and val not in p_before_set]
    
    # The optimal permutation starts with p_before, then m, then the rest.
    # The order of the rest doesn't matter.
    optimal_p = p_before + [m] + all_other_nums
    
    # Also add numbers from original input 'a' that were > initial_x
    p_large = sorted([val for val in a if val > initial_x])
    final_permutation = p_large + optimal_p

    # 6. Print the calculation steps
    x = initial_x
    print(f"Initial x = {x}")
    for val in final_permutation:
        new_x = x % val
        print(f"{x} % {val} = {new_x}")
        x = new_x
    print(f"Final x = {x}")


# Example usage:
# find_best_permutation([22, 35], 100)
# Expected output for [22, 35], 100 is 12 via permutation [22, 35]
# find_best_permutation([19, 20, 21], 100)
# Expected output for [19, 20, 21], 100 is 16 via permutation [21, 19, 20]
find_best_permutation([6, 7, 8], 27)
