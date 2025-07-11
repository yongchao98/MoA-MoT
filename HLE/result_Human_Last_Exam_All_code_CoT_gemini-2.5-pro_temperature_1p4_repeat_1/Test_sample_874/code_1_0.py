import math

# Memoization cache for the f_calc function to speed up computation
memo = {}

def get_canonical(a, b, c, d):
    """
    Finds the canonical representation of a tuple by checking all 8 cyclic
    permutations and reversals and picking the lexicographically smallest.
    This helps in effective memoization.
    """
    current = (a, b, c, d)
    # The initial tuple and its reversal are the first candidates for the minimum
    min_t = min(current, tuple(reversed(current)))

    # Iterate through the other 3 cyclic shifts
    for _ in range(3):
        current = (current[1], current[2], current[3], current[0])
        rev = tuple(reversed(current))
        min_t = min(min_t, current, rev)
    return min_t

def f_calc(a, b, c, d):
    """
    Calculates the length of the Ducci sequence for a given tuple (a, b, c, d).
    """
    # Use the canonical representation as the key for memoization
    key = get_canonical(a, b, c, d)
    if key in memo:
        return memo[key]

    # Base case for the recursion/termination
    if a == 0 and b == 0 and c == 0 and d == 0:
        return 1

    # The iterative process
    current = [a, b, c, d]
    count = 1
    # The sequence is guaranteed to terminate. The limit is a safeguard.
    limit = 100 
    for _ in range(limit):
        if all(x == 0 for x in current):
            break
        
        # Calculate the next square's values
        next_vals = [
            abs(current[0] - current[1]),
            abs(current[1] - current[2]),
            abs(current[2] - current[3]),
            abs(current[3] - current[0])
        ]
        current = next_vals
        count += 1
    
    # Store the result before returning
    memo[key] = count
    return count

def solve():
    """
    Main solver function to find the tuple and compute the final expression.
    """
    # 1. Generate the Tribonacci-like sequence with seed (1, 1, 2)
    U = [1, 1, 2]
    limit = 10_000_000
    while U[-1] <= limit:
        U.append(U[-1] + U[-2] + U[-3])
    
    # We can form tuples V_n as long as U[n+3] exists in our generated list
    max_n = len(U) - 4
    
    max_f = 0
    best_tuple = None
    min_sum = float('inf')

    # 2. Iterate through candidate tuples and find the optimal one
    for n in range(max_n + 1):
        a, b, c, d = U[n+3], U[n+2], U[n+1], U[n]
        
        # Calculate the length of the sequence
        length = f_calc(a, b, c, d)
        
        # Check if this tuple is better than the current best
        current_sum = a + b + c + d
        if length > max_f:
            max_f = length
            min_sum = current_sum
            best_tuple = (a, b, c, d)
        elif length == max_f:
            if current_sum < min_sum:
                min_sum = current_sum
                best_tuple = (a, b, c, d)

    a, b, c, d = best_tuple
    
    # 3. Compute the final result
    result = (a + b - c - d) % 1000
    
    print(f"The tuple (a,b,c,d) that maximizes f is {best_tuple}")
    print(f"The maximum value of f, M, is {max_f}")
    print(f"The sum a+b+c+d for this tuple is {min_sum}")
    print(f"The value of (a+b-c-d) is {a+b-c-d}")
    print(f"The final expression a + b - c - d mod 1000 evaluates to: {result}")
    
    # Present the final answer in the required format
    print(f"\nFinal Answer: {a} + {b} - {c} - {d} mod 1000 = {result}")
    print(f"<<<{result}>>>")

solve()
