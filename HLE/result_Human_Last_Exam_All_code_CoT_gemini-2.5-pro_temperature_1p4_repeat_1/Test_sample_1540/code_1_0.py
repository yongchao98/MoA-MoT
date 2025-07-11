import sys

def solve():
    """
    Reads two strings, compares them lexicographically ignoring case,
    and prints -1, 1, or 0.
    """
    try:
        s1 = sys.stdin.readline().strip().lower()
        s2 = sys.stdin.readline().strip().lower()

        if s1 < s2:
            print("-1")
        elif s1 > s2:
            print("1")
        else:
            print("0")
    except IndexError:
        # Handles empty input
        print("0")

solve()