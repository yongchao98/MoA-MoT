# A global counter to track the number of times the sort function is called.
call_count = 0

def swap(arr, i, j):
    """A helper function to swap two elements in a list."""
    arr[i], arr[j] = arr[j], arr[i]

def sort(a, n):
    """
    The provided sorting algorithm, instrumented to count recursive calls.
    """
    global call_count
    call_count += 1
    
    if n > 1:
        # First recursive call
        sort(a, n - 1)
        
        # This condition is always met for a reverse-sorted array,
        # leading to the worst-case behavior.
        if a[n - 1] < a[n - 2]:
            swap(a, n - 1, n - 2)
            # Second (expensive) recursive call
            sort(a, n - 1)

def demonstrate_complexity():
    """
    Runs the sort algorithm on worst-case inputs of increasing size
    to demonstrate its exponential time complexity.
    """
    print("Demonstrating the exponential growth of function calls for the worst-case scenario.")
    print("-" * 35)
    print(f"{'n':<5}{'Function Calls':<20}{'Expected (2^n - 1)':<20}")
    print("-" * 35)
    
    # We will test for n from 1 up to 16
    for n in range(1, 17):
        global call_count
        # Reset the counter for each run
        call_count = 0
        
        # The worst-case input is a reverse-sorted array
        worst_case_array = list(range(n, 0, -1))
        
        sort(worst_case_array, n)
        
        expected_calls = 2**n - 1
        
        print(f"{n:<5}{call_count:<20}{expected_calls:<20}")
        
demonstrate_complexity()