import sys

# The algorithm provided has a deep recursion, so we increase the limit.
# This is necessary for the demonstration, even for moderately small 'n'.
sys.setrecursionlimit(35000)

# We use a global variable to count the number of calls to the sort function.
# This will serve as our measure of the number of operations.
call_count = 0

def swap(arr, i, j):
    """A standard helper function to swap two elements in an array."""
    arr[i], arr[j] = arr[j], arr[i]

def sort_algorithm(a, n):
    """
    The sorting algorithm from the problem description, instrumented with a call counter.
    `a` is the list (array) and `n` is the size of the prefix to consider.
    """
    global call_count
    call_count += 1

    if n > 1:
        # 1. First recursive call on the prefix of size n-1
        sort_algorithm(a, n - 1)

        # 2. Compare the last two elements of the current prefix
        if a[n - 1] < a[n - 2]:
            swap(a, n - 1, a[n - 2])
            # 3. If a swap happened, make a second recursive call on the prefix of size n-1
            sort_algorithm(a, n - 1)

def analyze_complexity():
    """
    This function analyzes the time complexity by running the algorithm on worst-case
    inputs and explaining the resulting recurrence relation.
    """
    print("--- Analysis of the Sorting Algorithm's Worst-Case Time Complexity ---")
    print("\nThe worst-case input for this algorithm is a reverse-sorted array.")
    print("We will count the number of recursive calls for increasing array sizes (n).\n")
    print("Let's observe the growth pattern:")
    print("-" * 55)
    print(f"{'n':<5} | {'Recursive Calls':<20} | {'2^n':<20}")
    print("-" * 55)

    for n in range(1, 16):
        global call_count
        call_count = 0
        
        # Create a worst-case (reverse-sorted) array of size n
        worst_case_array = list(range(n, 0, -1))
        
        # Sort the array and count the calls
        sort_algorithm(worst_case_array, n)
        
        print(f"{n:<5} | {call_count:<20} | {2**n:<20}")

    print("-" * 55)

    print("\n--- Explanation of the Complexity ---")
    print("Let T(n) be the time complexity for an input of size n.")
    print("The function's execution in the worst case can be described by a recurrence relation.")
    print("\nIn the worst case (e.g., a reverse-sorted list), the `if` condition is always met,")
    print("triggering both recursive calls. The structure is:")
    print("  T(n) = T(n - 1)  (first recursive call)")
    print("       + C          (for the comparison and swap)")
    print("       + T(n - 1)  (second recursive call)")
    print("\nThis simplifies to the final equation for the recurrence relation:")
    print("  T(n) = 2 * T(n - 1) + C")

    print("\nSolving this relation:")
    print("  T(1) = C")
    print("  T(2) = 2*T(1) + C = 3C")
    print("  T(3) = 2*T(2) + C = 2*(3C) + C = 7C")
    print("  T(4) = 2*T(3) + C = 2*(7C) + C = 15C")
    print("The pattern is T(n) = (2^n - 1) * C.")
    print("\nThis means the number of operations grows exponentially with n.")
    print("Therefore, the worst-case time complexity is O(2^n).")

if __name__ == "__main__":
    analyze_complexity()