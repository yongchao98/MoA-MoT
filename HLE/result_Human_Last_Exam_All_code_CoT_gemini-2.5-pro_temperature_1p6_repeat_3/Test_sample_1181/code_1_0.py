def solve_solaris_riddle():
    """
    This function identifies and explains the answer to the Solaris trivia question.
    """
    choices = {
        'A': 'Kris',
        'B': 'Hari',
        'C': 'Snaut',
        'D': 'Sartorius',
        'E': 'Gibarian'
    }

    correct_option = 'A'
    correct_character = choices[correct_option]

    explanation = (
        f"The correct character is {correct_character}.\n"
        f"In the film, the protagonist Kris Kelvin expresses a deep sense of loss and shame for being disconnected from Earth.\n"
        "He specifically reflects on missing simple terrestrial experiences, symbolizing his profound homesickness and the psychological "
        "burden of his mission in the sterile environment of the space station."
    )
    
    print(explanation)

solve_solaris_riddle()