def compare_strings():
    """
    Reads two strings, compares them case-insensitively, and prints the result.
    -1 if the first is less than the second.
    1 if the first is greater than the second.
    0 if they are equal.
    """
    try:
        # Read the first string from input
        s1 = input()
        # Read the second string from input
        s2 = input()

        # Convert both strings to lowercase for case-insensitive comparison
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
        # Handle cases where input is not provided as expected
        pass

if __name__ == "__main__":
    compare_strings()
