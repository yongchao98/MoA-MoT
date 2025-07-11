# A global counter to track the number of calls
call_count = 0

def swap(arr, i, j):
    """Swaps two elements in an array."""
    arr[i], arr[j] = arr[j], arr[i]

def sort_instrumented(a, n):
    """The sorting algorithm with an instrumented call counter."""
    global call_count
    call_count += 1
    
    if n > 1:
        sort_instrumented(a, n - 1)
        if a[n - 1] < a[n - 2]:
            swap(a, n - 1, n - 2)
            sort_instrumented(a, n - 1)

def main():
    """
    Demonstrates the exponential complexity by running the algorithm
    on worst-case inputs of increasing size.
    """
    print("Analyzing the worst-case time complexity of the algorithm.")
    print("The recurrence relation is T(n) = 2 * T(n-1) + c.")
    print("This leads to O(2^n) complexity.")
    print("\n--- Demonstration by counting recursive calls ---")
    
    # Test for n from 1 to 10
    for n in range(1, 11):
        global call_count
        call_count = 0
        
        # Worst-case input: a reverse-sorted array
        arr = list(range(n, 0, -1))
        
        sort_instrumented(arr, n)
        
        # The number of calls for the recurrence C(n) = 2*C(n-1) + 1 is 2^n - 1.
        # Let's show this equation with the actual numbers.
        expected_calls = 2**n - 1
        
        print(f"\nFor n = {n}:")
        print(f"  Total calls: {call_count}")
        print(f"  Equation: 2^{n} - 1 = {expected_calls}")

if __name__ == "__main__":
    main()