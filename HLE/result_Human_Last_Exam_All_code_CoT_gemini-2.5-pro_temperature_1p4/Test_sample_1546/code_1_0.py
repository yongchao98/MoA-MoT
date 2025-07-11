import math

def get_s5_character_values():
    """
    Calculates the two irreducible characters of degree 4 for the symmetric group S5.
    """
    # The 7 conjugacy classes of S5 are defined by their cycle structure.
    # e.g., [2, 1, 1, 1] represents a transposition like (1 2).
    # The sum of elements in each sublist is 5.
    conjugacy_classes = [
        [1, 1, 1, 1, 1],  # identity
        [2, 1, 1, 1],      # transposition
        [2, 2, 1],          # double transposition
        [3, 1, 1],          # 3-cycle
        [3, 2],             # 3-cycle and 2-cycle
        [4, 1],             # 4-cycle
        [5],                # 5-cycle
    ]

    # These will store the character values for each of the 7 classes
    char_for_partition_4_1 = []
    char_for_partition_2_1_1_1 = []

    n = 5 # We are in S5

    for cycle_struct in conjugacy_classes:
        # --- Character 1: Corresponds to partition [4,1] ---
        # The character value is (number of fixed points) - 1.
        # Fixed points correspond to cycles of length 1.
        fixed_points = cycle_struct.count(1)
        val1 = fixed_points - 1
        char_for_partition_4_1.append(val1)

        # --- Character 2: Corresponds to partition [2,1,1,1] ---
        # This character is the product of the first character and the sign character.
        # The sign of a permutation is (-1)^(n - k), where k is the number of cycles.
        num_cycles = len(cycle_struct)
        sign = (-1)**(n - num_cycles)
        val2 = val1 * sign
        char_for_partition_2_1_1_1.append(val2)

    # The problem asks for the lists of numbers to be sorted in ascending order.
    char_for_partition_4_1.sort()
    char_for_partition_2_1_1_1.sort()

    # Print the two final lists, separated by a comma.
    print(char_for_partition_4_1, char_for_partition_2_1_1_1, sep=', ')

get_s5_character_values()