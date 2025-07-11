def solve():
    """
    Calculates the largest value k based on the formula 2 * (N - 2).
    N is the size of the square table, which is 2024.
    """
    N = 2024
    
    # The formula for the value k is 2 * (N - 2)
    k = 2 * (N - 2)

    # We print the steps of the calculation as requested.
    print("The calculation is based on the formula: k = 2 * (N - 2)")
    print(f"Given N = {N}:")
    print(f"k = 2 * ({N} - 2)")
    print(f"k = 2 * {N - 2}")
    print(f"k = {k}")

solve()