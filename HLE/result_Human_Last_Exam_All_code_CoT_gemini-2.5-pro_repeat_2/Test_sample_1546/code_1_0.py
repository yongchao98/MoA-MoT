def solve_s5_chars():
    """
    Calculates and prints the irreducible characters of degree 4 for the symmetric group S5.
    """
    # The 7 conjugacy classes of S5 are represented by the partitions of 5.
    n = 5
    partitions = [
        (1, 1, 1, 1, 1),  # Identity
        (2, 1, 1, 1),     # Transposition
        (2, 2, 1),       # Double transposition
        (3, 1, 1),         # 3-cycle
        (3, 2),            # (3-cycle)(2-cycle)
        (4, 1),            # 4-cycle
        (5,)               # 5-cycle
    ]

    # The first degree 4 character corresponds to the partition (4,1).
    # Its value is given by (number of fixed points) - 1.
    char1_values = []
    for p in partitions:
        # The number of fixed points is the number of 1s in the partition tuple.
        fixed_points = p.count(1)
        value = fixed_points - 1
        char1_values.append(value)

    # The second degree 4 character corresponds to the partition (2,1,1,1).
    # It is the product of the first character and the sign character.
    # The sign of a permutation is (-1)^(n - number of cycles).
    char2_values = []
    for i, p in enumerate(partitions):
        num_cycles = len(p)
        sign = (-1)**(n - num_cycles)
        value = char1_values[i] * sign
        char2_values.append(value)

    # Sort the character values in ascending order as requested.
    char1_values.sort()
    char2_values.sort()

    # Print the two lists, separated by a comma.
    # The output string must contain all numbers of the resulting lists.
    print(str(char1_values) + ", " + str(char2_values))

solve_s5_chars()