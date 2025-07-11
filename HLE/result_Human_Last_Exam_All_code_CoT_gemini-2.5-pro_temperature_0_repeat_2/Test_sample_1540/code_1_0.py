def solve():
    """
    Compares two strings lexicographically, ignoring case.
    """
    try:
        # Read the two strings from input
        s1 = input()
        s2 = input()

        # Convert both strings to the same case (e.g., lowercase)
        # for case-insensitive comparison.
        s1_lower = s1.lower()
        s2_lower = s2.lower()

        # Perform the lexicographical comparison
        if s1_lower < s2_lower:
            print("-1")
        elif s1_lower > s2_lower:
            print("1")
        else:
            print("0")

    except EOFError:
        # Handle case where there is no input
        pass

solve()