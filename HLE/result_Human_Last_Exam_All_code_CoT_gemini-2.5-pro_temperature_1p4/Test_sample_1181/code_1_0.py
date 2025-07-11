def solve_solaris_trivia():
    """
    This function identifies and prints the answer to the trivia question
    about the 1972 film 'Solaris'.
    """
    # The character who is ashamed to miss the sound of leaves rustling on Earth
    # is the protagonist, Kris Kelvin. This is a central theme of his character arc.
    answer_map = {
        "A": "Kris",
        "B": "Hari",
        "C": "Snaut",
        "D": "Sartorius",
        "E": "Gibarian"
    }

    correct_character = "Kris"
    correct_option = "A"

    print(f"In the 1972 Andrei Tarkovsky movie Solaris, the character ashamed to miss the sound of leaves rustling on Earth is {answer_map[correct_option]}.")
    print(f"This corresponds to answer choice {correct_option}.")

solve_solaris_trivia()