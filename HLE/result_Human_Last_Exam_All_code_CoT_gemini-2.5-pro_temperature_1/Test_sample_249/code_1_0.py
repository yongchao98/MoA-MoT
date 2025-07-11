import sys

def solve():
    """
    Calculates the minimum possible diameter of an undirected tree G
    with n+2 vertices and m leaves.
    """
    try:
        # Read n and m from command line arguments
        if len(sys.argv) != 3:
            print("Usage: python your_script.py <n> <m>")
            return

        n = int(sys.argv[1])
        m = int(sys.argv[2])

        if n <= 0 or m <= 0:
            print("Error: n and m must be positive integers.")
            return

        # A tree with more than 2 vertices must have at least 2 leaves.
        # If n+2 > 2 (i.e., n >= 1), we must have m >= 2.
        if n >= 1 and m < 2:
            print(f"Error: A tree with {n+2} vertices must have at least 2 leaves.")
            print("So, m must be >= 2 for n >= 1.")
            return

        # The number of leaves 'm' cannot exceed n+1 (for a star graph).
        if m > n + 1:
            print(f"Error: A tree with {n+2} vertices can have at most {n+1} leaves.")
            print("So, m cannot be greater than n+1.")
            return

        diameter = 0
        if m == n + 1:
            # Case 1: Star graph. This is possible for n>=1 (m>=2).
            diameter = 2
            equation_str = f"The minimum diameter is 2 since m = n + 1 (m={m}, n={n})."
        elif m == n:
            # Case 2: Bicentric graph. Diameter 3 is possible iff m=n.
            # This requires n>=2, as m=n and m>=2.
            if n < 2:
                 print(f"Error: The case m=n is only possible for n>=2 (since m>=2).")
                 return
            diameter = 3
            equation_str = f"The minimum diameter is 3 since m = n (m={m}, n={n})."
        else: # This covers m <= n - 1
            # Case 3: Caterpillar-like graph.
            diameter = n - m + 3
            equation_str = f"The minimum diameter is n - m + 3 = {n} - {m} + 3 = {diameter}."

        print(equation_str)

    except (ValueError, IndexError):
        print("Please provide two positive integers n and m as command-line arguments.")
        print("Usage: python your_script.py <n> <m>")

if __name__ == "__main__":
    solve()
