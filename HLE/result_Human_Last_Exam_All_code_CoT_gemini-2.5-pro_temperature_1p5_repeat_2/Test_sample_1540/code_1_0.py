def solve():
    """
    Compares two strings lexicographically, ignoring case.
    """
    try:
        # Read the first string from input
        first_string = input()

        # Read the second string from input
        second_string = input()

        # Convert both strings to lowercase for case-insensitive comparison
        s1_lower = first_string.lower()
        s2_lower = second_string.lower()

        # Compare the two lowercase strings and determine the output
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