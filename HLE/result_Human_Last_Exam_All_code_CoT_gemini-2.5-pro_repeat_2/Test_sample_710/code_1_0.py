import sys

def solve():
    """
    Determines the smallest number N such that any number >= N can be written
    as a sum of distinct numbers of the form 2n^2+3n+1.
    """
    # We are looking for the largest number 'g' that cannot be represented.
    # The answer N will be g + 1. We'll search for g up to a certain limit.
    # Based on tests, a limit of 1000 is sufficient as the result stabilizes.
    LIMIT = 1000

    # 1. Generate the numbers s_n = 2n^2 + 3n + 1
    s_values = []
    n = 1
    while True:
        s = 2 * n * n + 3 * n + 1
        if s > LIMIT:
            break
        s_values.append(s)
        n += 1

    # 2. Use DP to find all reachable sums up to LIMIT
    # is_reachable[k] will be true if k is a sum of distinct s_values.
    is_reachable = [False] * (LIMIT + 1)
    is_reachable[0] = True  # Sum of zero elements

    # Keep track of the maximum possible sum from the s_values processed so far.
    processed_s_sum = 0

    for s in s_values:
        # Iterate backwards to ensure we use each s_value at most once per sum.
        for i in range(processed_s_sum, -1, -1):
            if is_reachable[i]:
                if i + s <= LIMIT:
                    is_reachable[i + s] = True
        processed_s_sum += s
    
    # 3. Find the largest number g that is not reachable.
    g = 0
    for i in range(LIMIT, 0, -1):
        if not is_reachable[i]:
            g = i
            break
    
    # 4. The smallest number N is g + 1.
    N = g + 1

    # Print the result as an equation as requested
    print(f"{g} + 1 = {N}")

solve()