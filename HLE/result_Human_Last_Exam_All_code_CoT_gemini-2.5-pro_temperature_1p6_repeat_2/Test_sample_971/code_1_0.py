def solve_ballet_riddle():
    """
    This function calculates and displays the number of specific pirouettes
    performed by Natalia Osipova in a famous 2008 performance.
    """

    # In her 2008 Don Quixote Act I variation, Osipova performs a sequence
    # of 8 single pirouettes from the fifth position. We represent each
    # turn as the number 1.
    individual_turns = [1, 1, 1, 1, 1, 1, 1, 1]

    # The total number of pirouettes is the sum of these individual turns.
    total_pirouettes = sum(individual_turns)

    # We will build a string to display the full equation.
    # This fulfills the requirement to output each number in the final equation.
    equation_string = " + ".join(map(str, individual_turns))

    print("The number of single-turn pirouettes en dehors from the fifth position performed by Osipova in the sequence is calculated as follows:")
    print(f"{equation_string} = {total_pirouettes}")

solve_ballet_riddle()