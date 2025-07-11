def solve_solaris_question():
    """
    Analyzes character details from the 1972 movie Solaris to answer a specific question.
    """
    question = "In the 1972 Andrei Tarkovsky movie Solaris which character is ashamed to miss the sound of leaves rustling on Earth?"
    
    # Data structure containing information about the characters relevant to the question.
    characters = {
        'A': {
            'name': 'Kris Kelvin',
            'description': 'The protagonist who travels from a lush, natural Earth to the sterile Solaris station. In a key conversation, he reflects on his detachment from his home planet and expresses shame that he has nearly forgotten simple, natural sounds like the rustling of leaves.'
        },
        'B': {
            'name': 'Hari',
            'description': "A 'visitor' manifested from Kris's memories. As a creation of the Solaris ocean, she has no genuine memories of Earth to miss."
        },
        'C': {
            'name': 'Snaut',
            'description': "A cyberneticist on the station who is weary and philosophical but doesn't express this specific form of nostalgia."
        },
        'D': {
            'name': 'Sartorius',
            'description': 'A cold, scientifically-minded biologist who is focused on studying Solaris and shows little sentimentality for Earth.'
        },
        'E': {
            'name': 'Gibarian',
            'description': 'Kris\'s friend who is already deceased upon Kris\'s arrival and is only seen in a video recording.'
        }
    }

    correct_option = None
    reasoning = ""

    # Find the character who fits the description in the question.
    for option, details in characters.items():
        if "ashamed" in details['description'] and "rustling of leaves" in details['description']:
            correct_option = option
            reasoning = details['description']
            correct_character_name = details['name']
            break
            
    print(f"Question: {question}\n")
    print("Analyzing character motivations and dialogues...")
    print(f"The character who fits the description is {correct_character_name}.")
    print(f"\nReasoning: {reasoning}")
    print(f"\nTherefore, the correct option is {correct_option}.")

solve_solaris_question()