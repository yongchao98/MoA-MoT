import sys

# It's a deeply recursive function, so we might need to increase the recursion limit for larger N.
# sys.setrecursionlimit(2000)

# Global counter for operations
op_count = 0

def swap(arr, i, j):
    """Helper function to swap elements in an array."""
    arr[i], arr[j] = arr[j], arr[i]

def sort_and_count(a, n):
    """
    The sorting algorithm from the problem, modified to count comparisons.
    """
    global op_count
    if n > 1:
        sort_and_count(a, n - 1)
        
        op_count += 1 # Count the comparison as a basic operation
        if a[n - 1] < a[n - 2]:
            swap(a, n - 1, n - 2)
            sort_and_count(a, n - 1)

def demonstrate_complexity():
    """
    Demonstrates the algorithm's complexity by running it on worst-case inputs
    and showing that the operation count grows exponentially.
    """
    print("Analyzing the worst-case recurrence relation: T(n) = 2 * T(n-1) + 1")
    print("Assuming T(1) = 0 (no comparisons for n=1).")
    print("T(2) = 2 * T(1) + 1 = 2 * 0 + 1 = 1")
    print("T(3) = 2 * T(2) + 1 = 2 * 1 + 1 = 3")
    print("T(4) = 2 * T(3) + 1 = 2 * 3 + 1 = 7")
    print("T(5) = 2 * T(4) + 1 = 2 * 7 + 1 = 15")
    print("The number of comparisons follows the formula T(n) = 2^(n-1) - 1, which is O(2^n).\n")
    
    print("--- Verification with Code ---")
    # Test for n from 1 up to a reasonable limit to avoid long run times.
    for n in range(1, 13):
        global op_count
        op_count = 0
        
        # A reverse-sorted array is the worst-case input for this algorithm.
        worst_case_array = list(range(n, 0, -1))
        
        sort_and_count(worst_case_array, n)
        
        # The theoretical count of comparisons is 2^(n-1) - 1 for n > 1.
        theoretical_count = (2**(n - 1) - 1) if n > 1 else 0
        
        print(f"For n = {n}: counted {op_count} comparisons (Theoretical: {theoretical_count})")

if __name__ == '__main__':
    demonstrate_complexity()
