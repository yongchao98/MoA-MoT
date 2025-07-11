def main():
    """
    This script demonstrates the O(2^n) complexity of the given algorithm
    by showing how the number of recursive calls grows exponentially.
    """
    print("Analyzing the worst-case time complexity of the sorting algorithm.")
    print("The recurrence relation for the number of calls, C(n), is:")
    print("C(n) = 2 * C(n-1) + 1 for n > 1, with the base case C(1) = 1.")
    print("This relation solves to C(n) = 2^n - 1, which indicates O(2^n) complexity.\n")
    print("Let's calculate the number of calls for the first 10 values of n:")
    
    # We use the closed-form solution C(n) = 2^n - 1 for calculations
    for n in range(1, 11):
        if n == 1:
            print(f"For n = {n}: C(1) = 1 call")
        else:
            # C(n-1) = 2^(n-1) - 1
            prev_calls = 2**(n - 1) - 1
            # C(n) = 2^n - 1
            current_calls = 2**n - 1
            print(f"For n = {n}: C({n}) = 2 * C({n-1}) + 1 = 2 * {prev_calls} + 1 = {current_calls} calls")

    print("\nThe number of calls approximately doubles with each increment of n, a characteristic of exponential growth.")

if __name__ == "__main__":
    main()