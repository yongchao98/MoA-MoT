def solve_solaris_trivia():
    """
    This function identifies the correct character for the trivia question
    and prints the answer with an explanation.
    """
    question = "In the 1972 Andrei Tarkovsky movie Solaris which character is ashamed to miss the sound of leaves rustling on Earth?"
    
    answer_choices = {
        'A': 'Kris',
        'B': 'Hari',
        'C': 'Snaut',
        'D': 'Sartorius',
        'E': 'Gibarian'
    }

    correct_choice_letter = 'A'
    correct_choice_name = answer_choices[correct_choice_letter]

    explanation = (
        "In the film, the protagonist Kris Kelvin expresses this feeling. "
        "Overwhelmed by the sterile environment of the space station and the psychological phenomena "
        "he is experiencing, he feels a deep, nostalgic longing for Earth. He confesses that he is "
        "ashamed that he only now, when it's gone, truly misses a simple, natural sound like the rustling of leaves."
    )

    print(f"Question: {question}")
    print("-" * 20)
    print(f"The correct character is: {correct_choice_name}")
    print(f"Explanation: {explanation}")
    print("-" * 20)
    print(f"Final Answer Choice: {correct_choice_letter}")

solve_solaris_trivia()