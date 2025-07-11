import sys

def solve():
    """
    Solves the lexicographical string comparison problem.
    """
    
    # The memory calculation for the most memory-efficient C program is:
    # m = 101 (for char s1[101]) + 1 (for unsigned char i) + 4 (for int c2) = 106 bytes.
    # The numbers in the equation are 101, 1, and 4.

    # Read the two strings from standard input
    s1 = sys.stdin.readline().strip()
    s2 = sys.stdin.readline().strip()

    # Convert both strings to lowercase for case-insensitive comparison
    s1_lower = s1.lower()
    s2_lower = s2.lower()

    # Perform the comparison and print the result
    if s1_lower < s2_lower:
        print("-1")
    elif s1_lower > s2_lower:
        print("1")
    else:
        print("0")

solve()