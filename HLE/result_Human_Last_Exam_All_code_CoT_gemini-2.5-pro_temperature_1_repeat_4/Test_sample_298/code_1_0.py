import sys

# Increase recursion limit for the stirling2 function.
sys.setrecursionlimit(2000)

# Memoization cache for the Stirling numbers function.
memo = {}

def stirling2(n, k):
    """
    Calculates the Stirling number of the second kind, S(n, k),
    which is the number of ways to partition a set of n elements
    into k non-empty subsets.
    Uses recursion with memoization.
    """
    if (n, k) in memo:
        return memo[(n, k)]
    
    # Base cases
    if k < 0 or k > n:
        return 0
    if n == k:
        return 1
    if k == 0 and n > 0:
        return 0
    if k == 1 and n > 0:
        return 1

    # Recursive step
    result = k * stirling2(n - 1, k) + stirling2(n - 1, k - 1)
    memo[(n, k)] = result
    return result

def solve_cohomology(k):
    """
    Computes and prints the list of cohomology groups for M(k).
    """
    max_degree = 2 * (k - 1)
    cohomology_groups = []

    for d in range(max_degree + 1):
        if d % 2 != 0:
            # Odd degree cohomology groups are zero.
            cohomology_groups.append("0")
        else:
            j = d // 2
            # The rank is given by the Stirling number S(k, k-j).
            rank_k_param = k - j
            rank = stirling2(k, rank_k_param)

            # Format the group string based on the rank.
            if rank == 0:
                group_str = "0"
            elif rank == 1:
                group_str = "Z"
            else:
                group_str = f"Z^{rank}"
            cohomology_groups.append(group_str)
            
    # Print the final list in the specified format.
    print(f"[{', '.join(cohomology_groups)}]")

if __name__ == "__main__":
    # For the given problem, k = 7.
    k_val = 7
    solve_cohomology(k_val)