def solve_solaris_trivia():
    """
    This function solves the trivia question about the 1972 movie Solaris.
    """
    # The provided answer choices
    choices = {
        'A': 'Kris',
        'B': 'Hari',
        'C': 'Snaut',
        'D': 'Sartorius',
        'E': 'Gibarian'
    }

    # The correct answer is Kris Kelvin, the protagonist, who is overcome with
    # nostalgia for Earth, especially during the opening and closing scenes at his dacha.
    # His longing for simple, natural sounds like rustling leaves contrasts with the
    # sterile environment of the space station.
    correct_answer_key = 'A'
    correct_character_name = choices[correct_answer_key]

    print(f"In Andrei Tarkovsky's 1972 movie Solaris, the character who misses the sound of leaves rustling on Earth is {correct_answer_key}: {correct_character_name}.")

solve_solaris_trivia()