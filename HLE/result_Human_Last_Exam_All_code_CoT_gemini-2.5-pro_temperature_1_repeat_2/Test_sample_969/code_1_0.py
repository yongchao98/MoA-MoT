def solve_sequence():
    """
    This function solves the sequence puzzle by identifying a seed and reversing it.
    """
    # The given sequence
    sequence = [3, 2, 1, 2, 3, 3, 3, 2, 2]

    # The pattern appears to be an initial seed followed by repeating blocks.
    # The seed is the first four numbers, before the repeats start.
    seed = sequence[0:4]

    # The hypothesis is that the sequence continues by reversing the seed.
    next_four_elements = seed[::-1]

    # Print the full sequence including the original part and the next 4 elements.
    # The problem asks to output each number in the final equation.
    print(*sequence, *next_four_elements)

solve_sequence()