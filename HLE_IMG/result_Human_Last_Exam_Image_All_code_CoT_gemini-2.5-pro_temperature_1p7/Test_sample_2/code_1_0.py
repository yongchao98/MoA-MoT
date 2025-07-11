import sys

def analyze_sorting_algorithm():
    """
    Analyzes the provided recursive sorting algorithm by instrumenting it to count
    calls and running it on worst-case inputs.
    """
    
    # The recursion limit might be hit for n > ~20, so we can increase it if needed.
    # sys.setrecursionlimit(2000) 
    
    # A dictionary is used as a mutable counter object to be passed through recursion.
    counter = {'calls': 0}

    def sort_instrumented(a, n):
        """
        The provided sorting algorithm, instrumented to count recursive calls.
        """
        counter['calls'] += 1
        if n > 1:
            sort_instrumented(a, n - 1)
            # The comparison is a constant time operation.
            if a[n - 1] < a[n - 2]:
                # The swap is also a constant time operation.
                a[n - 1], a[n - 2] = a[n - 2], a[n - 1]
                sort_instrumented(a, n - 1)

    print("Analyzing the worst-case time complexity of the algorithm.")
    print("The recurrence relation for the number of calls, C(n), is:")
    
    # Outputting the numbers in the final equation as requested
    coefficient = 2
    n_minus_1_val = "n - 1"
    constant = 1 # for the call itself
    print(f"C(n) = {coefficient} * C({n_minus_1_val}) + {constant}")
    print("This recurrence solves to O(2^n). Let's verify experimentally.\n")

    max_n = 16
    for n in range(1, max_n + 1):
        # Reset counter for each run
        counter['calls'] = 0
        
        # A reverse-sorted array is a worst-case input for this algorithm.
        worst_case_array = list(range(n, 0, -1))
        
        # Run the instrumented sort
        sort_instrumented(worst_case_array, n)
        
        # The theoretical number of calls for C(n) = 2*C(n-1) + 1 with C(0)=0 is 2^n - 1.
        # But if we define C(1)=1, C(n) = 1 + 2*C(n-1) leads to C(n)=2^n-1 as well.
        theoretical_calls = 2**n - 1
        
        print(f"For n = {n:2d}: Actual calls = {counter['calls']:7d}, Theoretical calls (2^{n}-1) = {theoretical_calls:7d}")

analyze_sorting_algorithm()