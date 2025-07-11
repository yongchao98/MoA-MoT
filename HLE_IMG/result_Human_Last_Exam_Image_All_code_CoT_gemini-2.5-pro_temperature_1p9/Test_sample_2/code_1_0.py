def swap(arr, i, j):
    """Helper function to swap elements."""
    arr[i], arr[j] = arr[j], arr[i]

def sort_and_count(a, n, counter):
    """The given sort algorithm with a comparison counter."""
    # The base case is implicitly handled by the n > 1 check
    if n > 1:
        # First recursive call
        sort_and_count(a, n - 1, counter)

        # Increment counter for the comparison
        counter[0] += 1
        if a[n - 1] < a[n - 2]:
            swap(a, n - 1, n - 2)
            # Second recursive call (in the worst-case path)
            sort_and_count(a, n - 1, counter)

def analyze_complexity():
    """
    Analyzes the number of comparisons for the algorithm with worst-case inputs
    of increasing size n.
    """
    print("Analyzing the number of comparisons for the algorithm on worst-case inputs:")
    for n in range(1, 11):
        # A reverse-sorted array is a worst-case input.
        arr = list(range(n, 0, -1))
        # Use a list for the counter to make it mutable across function calls.
        counter = [0]
        sort_and_count(arr, n, counter)
        print(f"For n = {n}, number of comparisons is {counter[0]}")
    
    print("\nThe number of operations grows exponentially, which is characteristic of O(2^n).")
    print("For instance, the recurrence T(n) = 2*T(n-1) + c leads to exponential complexity.")

analyze_complexity()