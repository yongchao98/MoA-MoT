def solve_s5_chars():
    """
    This function calculates and prints the irreducible characters of degree 4
    for the symmetric group S5.
    """

    # The values of the first character of degree 4, corresponding to the partition [4,1].
    # These are derived from the formula (number of fixed points - 1) for each conjugacy class.
    char1_values = [4, 2, 1, 0, -1, 0, -1]

    # The values of the second character of degree 4, corresponding to the partition [2,1,1,1].
    # These are derived by multiplying the first character by the sign character.
    char2_values = [4, -2, 1, 0, -1, 0, 1]

    # Sort each list of character values in ascending order.
    sorted_char1 = sorted(char1_values)
    sorted_char2 = sorted(char2_values)

    # Print the two sorted lists, separated by a comma.
    print(f"{sorted_char1},{sorted_char2}")

solve_s5_chars()