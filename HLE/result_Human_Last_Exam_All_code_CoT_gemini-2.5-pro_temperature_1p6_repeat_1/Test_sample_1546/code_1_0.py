def solve():
    """
    Calculates and prints the irreducible characters of degree 4 for the symmetric group S5.
    """
    # The first character's values, corresponding to the partition (4,1).
    # These are derived from (number of fixed points) - 1 for each conjugacy class.
    char1_values = [4, 2, 1, 0, -1, 0, -1]

    # The second character's values, corresponding to the partition (2,1,1,1).
    # These are derived by multiplying the first character's values by the sign
    # of the corresponding conjugacy class.
    char2_values = [4, -2, 1, 0, -1, 1, 0]

    # Sort the lists of character values in ascending order, as requested.
    sorted_char1 = sorted(char1_values)
    sorted_char2 = sorted(char2_values)
    
    # Print the two lists, separated by a comma.
    # The format [list], [list] is achieved by converting lists to strings.
    print(str(sorted_char1) + ", " + str(sorted_char2))

solve()