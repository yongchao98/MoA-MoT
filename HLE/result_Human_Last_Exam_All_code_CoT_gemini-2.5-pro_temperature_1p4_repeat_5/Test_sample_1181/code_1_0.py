def solve_solaris_question():
    """
    This function identifies and prints the answer to a trivia question
    about the 1972 movie "Solaris".
    """
    question = "In the 1972 Andrei Tarkovsky movie Solaris which character is ashamed to miss the sound of leaves rustling on Earth?"

    answer_choices = {
        'A': 'Kris',
        'B': 'Hari',
        'C': 'Snaut',
        'D': 'Sartorius',
        'E': 'Gibarian'
    }

    correct_answer_key = 'E'
    
    print("Question:")
    print(question)
    print("\nAnswer Choices:")
    for key, value in answer_choices.items():
        print(f"{key}. {value}")
        
    # The character who says this is Dr. Gibarian, in the video message he leaves for Kris Kelvin.
    # In his final testament, he reflects on humanity's purpose in space and expresses shame for not having appreciated simple, earthly things more.
    
    correct_character_name = answer_choices[correct_answer_key]
    
    print(f"\nThe character ashamed to miss the sound of rustling leaves is {correct_answer_key}: {correct_character_name}.")

solve_solaris_question()