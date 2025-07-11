import sys

def solve_solaris_trivia():
    """
    This function identifies the character from the 1972 movie Solaris
    who is ashamed to miss the sound of leaves rustling on Earth.
    """
    # The character in question is the protagonist, the psychologist Kris Kelvin.
    # A central theme of the film is his profound homesickness and disconnection
    # from the natural world he left behind on Earth. His nostalgia is a powerful
    # force throughout the story, and the longing for simple sounds like
    # rustling leaves is a poignant expression of this.
    character_name = "Kris"
    
    # Mapping the character to the answer choices provided.
    answer_choices = {
        "A": "Kris",
        "B": "Hari",
        "C": "Snaut",
        "D": "Sartorius",
        "E": "Gibarian"
    }

    correct_choice = None
    for choice, name in answer_choices.items():
        if name == character_name:
            correct_choice = choice
            break

    print(f"The character from the 1972 movie 'Solaris' who is ashamed to miss the sound of leaves rustling on Earth is Kris Kelvin.")
    print(f"This corresponds to answer choice {correct_choice}.")

solve_solaris_trivia()