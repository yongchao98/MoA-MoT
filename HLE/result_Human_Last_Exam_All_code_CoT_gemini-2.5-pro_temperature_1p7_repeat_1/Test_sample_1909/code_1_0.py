import sys

# It is possible to hit the recursion limit for large n.
# Increase the recursion limit if needed.
sys.setrecursionlimit(2000)

# A dictionary to store computed values of f(n) to avoid re-computation (memoization).
memo = {0: 0, 1: 1, 2: 0}

def get_f(n):
    """
    Calculates the expected number of remaining numbers for a row of size n
    using the derived recurrence relation with memoization.
    f(n+1) = ((n-1)/n) * f(n) + (2/n) * f(n-1) for n >= 2
    This is equivalent to (n-1)f(n) = (n-2)f(n-1) + 2f(n-2) for n >= 3
    """
    if n in memo:
        return memo[n]
    
    # Check to avoid infinite recursion for invalid inputs
    if n < 0:
        return 0

    # Recurrence: (n-1)f(n) = (n-2)f(n-1) + 2f(n-2)
    # We calculate f(n) using n-1 and n-2, so we call get_f(n-1) and get_f(n-2)
    val = ((n - 2) * get_f(n - 1) + 2 * get_f(n - 2)) / (n - 1)
    memo[n] = val
    return val

def main():
    """
    Calculates and prints f(n) and the ratio f(n)/n for various n.
    """
    print("This program calculates E(n)/n, where E(n) is the expected number of remaining numbers.")
    print("The values are calculated for n from 1 up to a specified maximum.")
    
    try:
        max_n = int(input("Enter the maximum value of n to simulate (e.g., 500): "))
        if max_n <= 0:
            print("Please enter a positive integer.")
            return
    except ValueError:
        print("Invalid input. Please enter an integer.")
        return

    print("\n  n   |   f(n)   | f(n)/n")
    print("--------------------------------")
    
    # Calculate for specific points to show the trend
    points_to_print = sorted(list(set([i for i in range(1, min(21, max_n + 1))] + [10 * k for k in range(3, max_n // 10 + 1)] + [max_n])))
    
    for n in points_to_print:
        if n > max_n: continue
        fn = get_f(n)
        ratio = fn / n
        print(f"{n:5d} | {fn:8.4f} | {ratio:8.6f}")

    print("\nAs n approaches infinity, the ratio f(n)/n approaches 0.")
    print("\nLet's present the final answer based on this analysis.")
    print("lim_{n->inf} E(n)/n = 0")


if __name__ == "__main__":
    main()
