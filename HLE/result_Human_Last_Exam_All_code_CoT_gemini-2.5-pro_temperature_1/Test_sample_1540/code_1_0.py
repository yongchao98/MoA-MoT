import sys

def solve():
    """
    Reads two strings, compares them lexicographically ignoring case,
    and prints -1, 1, or 0.
    """
    # Read the first string from standard input and convert to lowercase.
    s1 = sys.stdin.readline().strip().lower()

    # Read the second string from standard input and convert to lowercase.
    s2 = sys.stdin.readline().strip().lower()

    # Perform the comparison.
    if s1 < s2:
        print("-1")
    elif s1 > s2:
        print("1")
    else:
        print("0")

solve()