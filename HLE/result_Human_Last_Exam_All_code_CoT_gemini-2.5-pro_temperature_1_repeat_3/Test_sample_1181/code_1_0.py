def solve_solaris_question():
    """
    This script solves a trivia question about the 1972 film 'Solaris'.
    It identifies the character who laments missing the sound of rustling leaves.
    """

    # The multiple-choice options provided in the question
    answer_choices = {
        'A': 'Kris',
        'B': 'Hari',
        'C': 'Snaut',
        'D': 'Sartorius',
        'E': 'Gibarian'
    }

    # The correct character is Dr. Snaut. In a conversation with Kris Kelvin,
    # he cynically reflects on their mission and their disconnection from Earth.
    # He expresses a longing for simple, human, terrestrial experiences,
    # specifically mentioning being ashamed to miss the sound of rustling leaves,
    # which symbolizes their profound isolation and loss of humanity.
    correct_answer_letter = 'C'
    correct_character_name = answer_choices[correct_answer_letter]

    # Print the explanation
    print("In Andrei Tarkovsky's 1972 movie 'Solaris', the character who expresses this sentiment is Dr. Snaut.")
    print(f"In a moment of philosophical reflection, {correct_character_name} tells the protagonist, Kris, about the things they've lost by being in space.")
    print("He is ashamed to have forgotten such a simple, natural sound from Earth, highlighting the psychological toll of their mission and their separation from their home world.")
    print("\nTherefore, the correct answer is C.")

solve_solaris_question()