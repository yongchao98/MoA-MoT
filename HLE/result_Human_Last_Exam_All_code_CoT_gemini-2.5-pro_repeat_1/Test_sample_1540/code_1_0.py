import sys

def solve():
    """
    This function reads two strings, compares them lexicographically ignoring case,
    and prints -1, 1, or 0 as a result.
    """
    # Read the first string from standard input.
    try:
        s1 = sys.stdin.readline().strip()
    except (IOError, IndexError):
        s1 = ""

    # Read the second string from standard input.
    try:
        s2 = sys.stdin.readline().strip()
    except (IOError, IndexError):
        s2 = ""

    # Convert both strings to lowercase for case-insensitive comparison.
    s1_lower = s1.lower()
    s2_lower = s2.lower()

    # Perform the lexicographical comparison and determine the result.
    if s1_lower < s2_lower:
        result = -1
    elif s1_lower > s2_lower:
        result = 1
    else:
        result = 0

    # The problem statement mentions an equation, which is not applicable here.
    # We will just print the final result as requested.
    # Final Result: -1
    if result == -1:
        print("-1")
    # Final Result: 1
    elif result == 1:
        print("1")
    # Final Result: 0
    else:
        print("0")

solve()