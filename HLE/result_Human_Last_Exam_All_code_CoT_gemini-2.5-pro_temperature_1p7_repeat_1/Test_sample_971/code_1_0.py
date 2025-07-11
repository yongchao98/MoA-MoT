def solve_ballet_question():
    """
    This function calculates the number of turns in a specific sequence
    of Kitri's Act I variation in Don Quixote.

    The famous diagonal turning sequence in this variation is traditionally
    composed of 8 turns. We can represent this as a calculation.
    """

    # The sequence can be broken down into two symmetrical phrases.
    number_of_phrases = 2

    # Each phrase contains a set of 4 turns.
    turns_per_phrase = 4

    # Calculate the total number of turns.
    total_turns = number_of_phrases * turns_per_phrase

    # Output the final equation and the result.
    print("In the iconic diagonal turning sequence of the Act I variation, the number of turns is calculated as follows:")
    print(f"{number_of_phrases} * {turns_per_phrase} = {total_turns}")

solve_ballet_question()