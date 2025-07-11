import sys

def solve():
    """
    Solves the Petya and Strings problem.
    """
    try:
        # Read the two strings from standard input
        s1 = sys.stdin.readline().strip()
        s2 = sys.stdin.readline().strip()

        # Convert both strings to lowercase for case-insensitive comparison
        s1_lower = s1.lower()
        s2_lower = s2.lower()

        # Compare the two lowercase strings lexicographically
        if s1_lower < s2_lower:
            print("-1")
        elif s1_lower > s2_lower:
            print("1")
        else:
            print("0")

    except Exception as e:
        # Handle potential errors, though problem constraints suggest valid input
        print(f"An error occurred: {e}", file=sys.stderr)

solve()
