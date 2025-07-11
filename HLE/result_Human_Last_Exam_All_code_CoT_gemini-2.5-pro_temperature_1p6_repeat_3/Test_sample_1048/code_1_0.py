import sys

# It's a recursive problem, so we might hit the recursion limit.
sys.setrecursionlimit(2000)

# Memoization table to store results of subproblems
# The key will be a tuple (x, mask), and the value will be (max_value, best_permutation)
memo = {}
initial_x = 0
num_list = []
M = 0

def solve(x, mask):
    """
    Recursively finds the maximum final x and the permutation that achieves it.
    
    Args:
        x (int): The current value.
        mask (int): A bitmask representing the set of available numbers from num_list.
    
    Returns:
        tuple: A pair containing (maximum final value, corresponding permutation).
    """
    if mask == 0:
        # Base case: no more numbers to apply
        return (x, [])

    # We only memoize states where x is smaller than the maximum number in the list.
    # The initial x can be very large, and memoizing it is inefficient.
    is_memoizable = x < M
    if is_memoizable and (x, mask) in memo:
        return memo[(x, mask)]

    best_val = -1
    best_perm = []

    for i in range(len(num_list)):
        # Check if the i-th number is in the current available set (mask)
        if (mask >> i) & 1:
            # Recursively call solve with the new value of x and updated mask
            # new_mask removes the i-th number
            new_mask = mask ^ (1 << i)
            val, perm = solve(x % num_list[i], new_mask)

            if val > best_val:
                best_val = val
                # The new permutation is the current number plus the permutation from the subproblem
                best_perm = [num_list[i]] + perm
    
    if is_memoizable:
        memo[(x, mask)] = (best_val, best_perm)
        
    return (best_val, best_perm)

def find_best_permutation(x, a):
    """
    Main function to kick off the calculation.
    """
    global initial_x, num_list, M, memo
    initial_x = x
    # Remove duplicates as they don't add new paths (x % k % k is same as x % k)
    num_list = sorted(list(set(a))) 
    if not num_list:
        print(f"{initial_x} = {initial_x}")
        return

    M = max(num_list) if num_list else 0
    memo.clear()

    full_mask = (1 << len(num_list)) - 1
    
    final_x, permutation = solve(initial_x, full_mask)
    
    equation_str = str(initial_x)
    for num in permutation:
        equation_str += f" % {num}"
    
    equation_str += f" = {final_x}"
    
    print(equation_str)

# Example Usage:
# Given a list of n numbers and another number x
n_numbers = [10, 20, 3]
x_val = 25
find_best_permutation(x_val, n_numbers)

n_numbers_2 = [300, 1000]
x_val_2 = 1400
find_best_permutation(x_val_2, n_numbers_2)

n_numbers_3 = [2, 3, 100]
x_val_3 = 123
find_best_permutation(x_val_3, n_numbers_3)
