import sys

# It's possible to reach python's recursion limit for deep sequences.
# Increase the recursion limit for safety, although memoization should help.
sys.setrecursionlimit(2000)

# Memoization cache for the f() function
memo_f = {}

def f(a, b, c, d):
    """
    Calculates the length of the Ducci sequence for a given tuple (a, b, c, d).
    Uses memoization to store and retrieve results for previously seen tuples.
    """
    # Create a tuple to use as a dictionary key
    start_tuple = (a, b, c, d)
    
    # Return cached result if available
    if start_tuple in memo_f:
        return memo_f[start_tuple]

    # Base case: the all-zero tuple has a sequence length of 1
    if start_tuple == (0, 0, 0, 0):
        return 1

    # Calculate the next tuple in the sequence
    next_tuple = (abs(a - b), abs(b - c), abs(c - d), abs(d - a))
    
    # The length of the current sequence is 1 + the length of the sequence from the next tuple.
    # Recursively call f() for the next tuple and add 1.
    result = 1 + f(*next_tuple)
    
    # Cache the result before returning
    memo_f[start_tuple] = result
    return result

def solve():
    """
    Finds the optimal tuple and computes the final expression.
    """
    # --- Step 1: Generate Tribonacci numbers up to the limit ---
    limit = 10_000_000
    trib_nums = [0, 1, 1]
    while True:
        next_trib = trib_nums[-1] + trib_nums[-2] + trib_nums[-3]
        if next_trib > limit:
            break
        trib_nums.append(next_trib)
    
    # --- Step 2: Find the 'n' that maximizes f for the Tribonacci family tuple ---
    max_len = 0
    best_n = -1 # The index in trib_nums that gives the max length

    # We iterate through possible n values. T_n is trib_nums[n].
    # We need a tuple (T_n, T_{n-1}, T_{n-2}, 0), so n must be at least 2.
    for n in range(2, len(trib_nums)):
        a = trib_nums[n]
        b = trib_nums[n-1]
        c = trib_nums[n-2]
        d = 0
        
        length = f(a, b, c, d)
        
        # If we find a new maximum length, update max_len and best_n.
        # We want the smallest n for the minimal sum, so we don't update if length == max_len.
        if length > max_len:
            max_len = length
            best_n = n

    # --- Step 3: Determine the optimal tuple and calculate the result ---
    a = trib_nums[best_n]
    b = trib_nums[best_n - 1]
    c = trib_nums[best_n - 2]
    d = 0
    
    print(f"The maximum length M found is: {max_len}")
    print(f"This is achieved with n = {best_n}, giving the tuple with the smallest sum.")
    print(f"a = {a}")
    print(f"b = {b}")
    print(f"c = {c}")
    print(f"d = {d}")
    
    result = a + b - c - d
    final_answer = result % 1000
    
    print(f"\nThe expression to compute is (a + b - c - d) mod 1000")
    print(f"Calculation: ({a} + {b} - {c} - {d}) mod 1000")
    print(f"= ({a + b} - {c}) mod 1000")
    print(f"= {result} mod 1000")
    print(f"= {final_answer}")
    
    print(f"\nFinal Answer: {final_answer}")
    return final_answer

# Execute the solver and print the final answer in the desired format
final_result = solve()
print(f"<<<{final_result}>>>")
