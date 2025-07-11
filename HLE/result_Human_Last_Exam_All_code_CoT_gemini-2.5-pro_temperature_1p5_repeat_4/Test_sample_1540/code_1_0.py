import sys

def solve():
    """
    Reads two strings, compares them lexicographically ignoring case,
    and prints -1, 1, or 0.
    """
    # Read the first string from standard input.
    s1 = sys.stdin.readline().strip()
    
    # Read the second string from standard input.
    s2 = sys.stdin.readline().strip()
    
    # Convert both strings to lowercase for case-insensitive comparison.
    s1_lower = s1.lower()
    s2_lower = s2.lower()
    
    # Perform the lexicographical comparison and print the result.
    if s1_lower < s2_lower:
        print("-1")
    elif s1_lower > s2_lower:
        print("1")
    else:
        print("0")

solve()