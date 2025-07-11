import sys

# It's necessary to increase the recursion limit for demonstrating the rapid growth.
sys.setrecursionlimit(2000)

# Global variable to count the number of calls to the sort function.
call_count = 0

def swap(arr, i, j):
    """Helper function to swap two elements in a list."""
    arr[i], arr[j] = arr[j], arr[i]

def sort(a, n):
    """The sorting algorithm from the problem, with a call counter."""
    global call_count
    call_count += 1
    
    if n > 1:
        sort(a, n - 1)
        if a[n - 1] < a[n - 2]:
            swap(a, n - 1, n - 2)
            sort(a, n - 1)

def analyze_complexity():
    """
    Analyzes and demonstrates the time complexity of the sort function.
    """
    print("Analyzing the worst-case time complexity of the given sorting algorithm.")
    print("The worst case occurs with a reverse-sorted list.\n")
    print("The recurrence relation for the number of operations T(n) in the worst case is:")
    print("T(n) = T(n-1) + T(n-1) + O(1)")
    print("T(n) = 2 * T(n-1) + c\n")
    print("This recurrence resolves to O(2^n). Let's verify this by counting function calls.")
    print("-" * 60)
    print(f"{'n':<5} | {'Actual Calls':<20} | {'Theoretical (2^n - 1)':<25}")
    print("-" * 60)
    
    for n in range(1, 16):
        global call_count
        call_count = 0
        
        # A reverse-sorted list is the worst case for this algorithm.
        worst_case_array = list(range(n, 0, -1))
        
        # Call the sort function.
        sort(worst_case_array, n)
        
        # The number of calls follows the pattern T(n) = 2T(n-1) + 1, with T(1)=1.
        # This resolves to T(n) = 2^n - 1.
        theoretical_calls = 2**n - 1
        
        print(f"{n:<5} | {call_count:<20} | 2^{n} - 1 = {theoretical_calls}")

if __name__ == "__main__":
    analyze_complexity()
