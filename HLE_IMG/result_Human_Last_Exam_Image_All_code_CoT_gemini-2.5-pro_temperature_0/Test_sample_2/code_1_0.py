import sys

# Use a list to pass the counter by reference, allowing it to be modified
# by the recursive function calls.
call_count = [0]

def swap(arr, i, j):
    """A helper function to swap two elements in an array."""
    arr[i], arr[j] = arr[j], arr[i]

def sort(a, n):
    """
    The recursive sorting algorithm from the problem.
    This implementation also counts the number of times it is called.
    """
    # Increment the call counter each time the function is entered.
    call_count[0] += 1

    if n > 1:
        # First recursive call on the first n-1 elements.
        sort(a, n - 1)
        
        # Compare the last two elements of the current sub-array.
        # Note: Python uses 0-based indexing, so a[n-1] is the nth element.
        if a[n - 1] < a[n - 2]:
            # Swap if the last element is smaller than the one before it.
            swap(a, n - 1, n - 2)
            # Second recursive call, which leads to exponential complexity.
            sort(a, n - 1)

def main():
    """
    Main function to demonstrate the algorithm's complexity by counting calls.
    """
    # Increase recursion limit for larger n, though we'll keep n small.
    # sys.setrecursionlimit(3000)
    
    print("Analyzing the time complexity of the given recursive algorithm.")
    print("The recurrence relation for the number of calls in the worst case is T(n) = 2*T(n-1) + 1.")
    print("This solves to T(n) = 2^n - 1, which is O(2^n).")
    print("Let's verify this by counting the calls for different values of n.")
    print("-" * 60)
    print(f"{'n':<5}{'Actual Calls':<20}{'Expected (2^n - 1)':<25}")
    print("-" * 60)

    # Test for n from 1 to 15.
    for n in range(1, 16):
        # Reset the counter for each run.
        global call_count
        call_count = [0]
        
        # A reverse-sorted array is a worst-case input for this algorithm.
        arr = list(range(n, 0, -1))
        
        # Run the sort. We are interested in the call count, not the sorted result.
        sort(arr, n)
        
        # The theoretical number of calls is 2^n - 1.
        expected_calls = 2**n - 1
        
        print(f"{n:<5}{call_count[0]:<20}{expected_calls:<25}")

if __name__ == "__main__":
    main()