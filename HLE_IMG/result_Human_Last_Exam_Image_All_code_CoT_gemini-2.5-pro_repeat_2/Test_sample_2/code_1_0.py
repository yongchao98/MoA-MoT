import sys

# It's necessary to increase the recursion limit for this inefficient algorithm,
# even for small n.
sys.setrecursionlimit(3000)

# Global counter for function calls to measure complexity
call_count = 0

def swap(arr, i, j):
    """Swaps two elements in a list."""
    arr[i], arr[j] = arr[j], arr[i]

def sort_instrumented(a, n):
    """
    The sorting algorithm from the problem, instrumented with a counter
    to track the number of recursive calls.
    """
    global call_count
    call_count += 1
    
    if n > 1:
        sort_instrumented(a, n - 1)
        # In Python, list indices are 0-based.
        # The pseudocode's a[n-1] and a[n-2] correspond to the same indices.
        if a[n - 1] < a[n - 2]:
            swap(a, n - 1, n - 2)
            sort_instrumented(a, n - 1)

def analyze_complexity():
    """
    Analyzes the complexity by running the sort on worst-case inputs
    and printing the number of recursive calls.
    """
    print("### Empirical Verification ###")
    print("\nRunning the algorithm on a worst-case input (a reverse-sorted array)")
    print("and counting the number of recursive calls.")
    print("-" * 50)
    print(f"{'n':<5}{'Recursive Calls':<20}{'2^n - 1':<20}")
    print("-" * 50)
    
    # We can't test large n due to exponential complexity and recursion depth
    for n in range(1, 16):
        global call_count
        call_count = 0
        
        # Create a worst-case array (e.g., for n=5, arr = [4, 3, 2, 1, 0])
        worst_case_array = list(range(n - 1, -1, -1))
        
        try:
            # We pass the list and its length, as per the pseudocode
            sort_instrumented(worst_case_array, n)
            expected_calls = 2**n - 1
            print(f"{n:<5}{call_count:<20}{expected_calls:<20}")
            
            # Sanity check to ensure the array is actually sorted
            if worst_case_array != list(range(n)):
                print(f"Error: Array not sorted for n={n}")

        except RecursionError:
            print(f"{n:<5}{'Recursion limit reached':<20}")
            break

def main():
    """Main function to explain the analysis and run the verification."""
    print("### Analysis of the Algorithm's Time Complexity ###")
    print("\nTo find the time complexity, we derive and solve its recurrence relation.")
    
    print("\nStep 1: Formulate the Recurrence Relation")
    print("Let T(n) be the number of operations for an input of size n in the worst case.")
    print("The worst case occurs when the condition `a[n - 1] < a[n - 2]` is always true,")
    print("triggering two recursive calls for a problem of size n-1.")
    
    print("\nThe recurrence relation is: T(n) = T(n-1) + T(n-1) + c")
    print("T(n) = 2 * T(n-1) + c  (where c is a constant for the compare/swap operations)")
    print("The base case is T(1) = c'.")
    
    print("\nStep 2: Solve the Recurrence Relation")
    print("We can solve this by expansion, showing each number in the equation:")
    print("For T(1) = 1, c = 1 (counting calls):")
    print("T(1) = 1                  = 2^1 - 1")
    print("T(2) = 2 * T(1) + 1 = 2 * 1 + 1 = 3  = 2^2 - 1")
    print("T(3) = 2 * T(2) + 1 = 2 * 3 + 1 = 7  = 2^3 - 1")
    print("T(4) = 2 * T(3) + 1 = 2 * 7 + 1 = 15 = 2^4 - 1")
    print("\nThis pattern shows that T(n) = 2^n - 1.")
    print("Since the number of operations is proportional to 2^n, the time complexity is O(2^n).")
    
    # Run the code to provide empirical evidence
    analyze_complexity()

if __name__ == "__main__":
    main()