def solve_solaris_question():
    """
    Analyzes the characters of the 1972 film "Solaris" to answer a specific trivia question.
    """
    characters = {
        'Kris': {
            'role': "The protagonist, a psychologist sent to the space station.",
            'connection_to_earth': "Deeply melancholic and nostalgic for Earth. His journey is defined by his memories, guilt, and longing for his home planet and his past life. He explicitly laments having forgotten simple earthly experiences.",
            'is_correct': True
        },
        'Hari': {
            'role': "A 'visitor', a physical manifestation of Kris's deceased wife.",
            'connection_to_earth': "She is not human and has no memories of Earth beyond what she derives from Kris's mind. Therefore, she cannot miss something she never truly experienced.",
            'is_correct': False
        },
        'Snaut': {
            'role': "A cynical scientist on the station.",
            'connection_to_earth': "He has become detached and philosophical after years in space. While weary, he does not express the specific, shameful nostalgia mentioned in the question.",
            'is_correct': False
        },
        'Sartorius': {
            'role': "A cold, rational scientist on the station.",
            'connection_to_earth': "He is focused entirely on the scientific problem of Solaris and views emotions and nostalgia as hindrances. He shows no longing for Earth.",
            'is_correct': False
        },
        'Gibarian': {
            'role': "Kris's colleague who died before his arrival.",
            'connection_to_earth': "Appears in a video message expressing distress, but the specific line about shame and missing the sound of leaves is spoken by Kris.",
            'is_correct': False
        }
    }

    answer = None
    reasoning = ""
    
    print("Finding the character who missed the sound of rustling leaves...\n")

    for name, data in characters.items():
        if data['is_correct']:
            answer = name
            reasoning = data['connection_to_earth']
            break

    if answer:
        print(f"The correct character is Kris.")
        print("\nReasoning:")
        print(reasoning)
    else:
        print("Could not determine the correct answer.")

solve_solaris_question()