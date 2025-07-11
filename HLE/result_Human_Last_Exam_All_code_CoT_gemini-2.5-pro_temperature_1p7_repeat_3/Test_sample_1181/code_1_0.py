def find_solaris_character():
    """
    This function identifies the character from the 1972 film 'Solaris'
    who expressed shame at missing the sound of rustling leaves.
    """
    
    # A dictionary of the characters provided in the answer choices.
    character_options = {
        "A": "Kris",
        "B": "Hari",
        "C": "Snaut",
        "D": "Sartorius",
        "E": "Gibarian"
    }
    
    # Dr. Snaut is the character who delivers a philosophical monologue
    # about humanity's homesickness and the need for a "mirror" rather than new worlds.
    # It is he who is ashamed to admit missing simple, earthly things like the rustling of leaves.
    correct_answer_key = "C"
    correct_character_name = character_options[correct_answer_key]
    
    print("In Andrei Tarkovsky's 1972 movie Solaris, the character who is ashamed to miss the sound of leaves rustling on Earth is:")
    print(f"Answer: {correct_answer_key}. {correct_character_name}")

find_solaris_character()