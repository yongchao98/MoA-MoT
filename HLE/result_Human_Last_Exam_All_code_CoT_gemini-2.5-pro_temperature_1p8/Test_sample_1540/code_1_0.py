import sys

def solve():
    """
    Solves the Petya and Strings problem.
    """
    try:
        # Read the first string from input
        s1 = sys.stdin.readline().strip()
        # Read the second string from input
        s2 = sys.stdin.readline().strip()

        # Convert both strings to lowercase for case-insensitive comparison
        s1_lower = s1.lower()
        s2_lower = s2.lower()

        # Perform lexicographical comparison
        if s1_lower < s2_lower:
            print("-1")
        elif s1_lower > s2_lower:
            print("1")
        else:
            print("0")

    except IndexError:
        # Handle cases with fewer than two lines of input gracefully
        pass

solve()