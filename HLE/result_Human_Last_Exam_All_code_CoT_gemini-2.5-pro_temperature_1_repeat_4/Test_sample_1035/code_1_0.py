import math

memo = {}

def get_cost_and_split(n):
    """
    Calculates the minimum cost C(n) and the optimal split point using dynamic programming.
    Returns a tuple (cost, split_point).
    """
    if n in memo:
        return memo[n]

    if n <= 1:
        return (0, 0)

    # Check if n is a power of 2
    if (n > 0) and (n & (n - 1) == 0):
        k = int(math.log2(n))
        cost = k * n
        # For powers of 2, the canonical split is n/2.
        split = n // 2
        memo[n] = (cost, split)
        return cost, split

    min_cost = float('inf')
    best_split = -1

    for m in range(1, n // 2 + 1):
        cost_m, _ = get_cost_and_split(m)
        cost_nm, _ = get_cost_and_split(n - m)
        current_cost = n + cost_m + cost_nm
        
        if current_cost < min_cost:
            min_cost = current_cost
            best_split = m

    memo[n] = (min_cost, best_split)
    return min_cost, best_split

def solve():
    """
    Solves the problem for n=35 and prints the result in the required format.
    """
    n = 35
    
    # Populate the memoization table up to n
    for i in range(1, n + 1):
        get_cost_and_split(i)
        
    total_cost, split_m = get_cost_and_split(n)
    split_nm = n - split_m
    
    cost_m, _ = get_cost_and_split(split_m)
    cost_nm, _ = get_cost_and_split(split_nm)
    
    print(f"The minimum number of vertices is found by solving the recurrence C(n) = n + C(m) + C(n-m).")
    print(f"For n = {n}, the optimal split is m = {split_m} and n-m = {split_nm}.")
    print(f"The costs for the subproblems are C({split_m}) = {cost_m} and C({split_nm}) = {cost_nm}.")
    print(f"The final calculation is:")
    print(f"{n} + {cost_m} + {cost_nm} = {total_cost}")

solve()