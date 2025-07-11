def find_solaris_character():
    """
    This function identifies and prints the character from the 1972 movie "Solaris"
    who is ashamed to miss the sound of rustling leaves on Earth.
    """
    # The list of characters from the movie provided in the options.
    characters = {
        "A": "Kris",
        "B": "Hari",
        "C": "Snaut",
        "D": "Sartorius",
        "E": "Gibarian"
    }
    
    # The correct character is Dr. Snaut.
    correct_option = "C"
    character_name = characters[correct_option]
    
    # Print the explanation and the answer.
    print(f"In the 1972 film 'Solaris', the character who expresses shame for missing earthly sounds like rustling leaves is Dr. {character_name}.")
    print(f"The correct answer choice is: {correct_option}")

find_solaris_character()