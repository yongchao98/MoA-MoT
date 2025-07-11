# Global counter to track the number of recursive calls
call_counter = 0

def swap(arr, i, j):
    """Helper function to swap two elements in a list."""
    arr[i], arr[j] = arr[j], arr[i]

def sort(a, n):
    """
    The sorting algorithm from the problem, instrumented with a call counter.
    """
    global call_counter
    call_counter += 1

    if n > 1:
        # First recursive call
        sort(a, n - 1)
        
        # In Python, arrays are 0-indexed, so we access elements
        # at indices n-1 and n-2.
        if a[n - 1] < a[n - 2]:
            swap(a, n - 1, n - 2)
            # Second recursive call (worst-case path)
            sort(a, n - 1)

def analyze_complexity():
    """
    Runs the algorithm on worst-case inputs and prints the number of calls
    to demonstrate the O(2^n) complexity.
    """
    print("--- Analysis of the Algorithm's Time Complexity ---")
    print("The worst-case recurrence relation is T(n) = 2 * T(n-1) + c.")
    print("This results in an exponential time complexity of O(2^n).")
    print("\nDemonstration using a worst-case input (reverse-sorted array):")
    print("We will count the total number of calls to the sort function.")
    
    print("-" * 55)
    print(f"{'Input Size (n)':<15} | {'Theoretical Calls (2^n - 1)':<25} | {'Actual Calls':<15}")
    print("-" * 55)
    
    global call_counter
    for n in range(1, 16):
        # Reset the counter for each run
        call_counter = 0
        
        # A reverse-sorted list is a worst-case input
        worst_case_array = list(range(n, 0, -1))
        
        # Run the algorithm
        sort(worst_case_array, n)
        
        # The total number of calls in the worst case is 2^n - 1
        expected_calls = 2**n - 1
        
        print(f"{n:<15} | {expected_calls:<25} | {call_counter:<15}")

# Run the analysis
analyze_complexity()
