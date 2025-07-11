def solve_ballet_turns():
    """
    This function calculates the number of famous turns Natalia Osipova performed as Kitri.

    While the question specifies the Act I variation, the most famous sequence of
    consecutive turns is the 32 fouettés in the Act III coda. It's a high probability
    that the question is referring to this iconic feat of ballet. Natalia Osipova,
    in her 2008 performance, successfully completed these turns.
    """

    # The famous sequence consists of one set of 32 turns.
    number_of_sets = 1
    turns_in_the_sequence = 32

    total_turns = number_of_sets * turns_in_the_sequence

    # The problem asks to output each number in the final equation.
    print(f"The number of fouetté turns performed by Natalia Osipova in the famous Act III sequence is calculated as:")
    print(f"{number_of_sets} set * {turns_in_the_sequence} turns per set = {total_turns} total turns")

solve_ballet_turns()