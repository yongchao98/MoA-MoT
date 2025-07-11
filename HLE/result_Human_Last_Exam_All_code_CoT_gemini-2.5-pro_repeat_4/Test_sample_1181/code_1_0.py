import sys

def solve_solaris_question():
    """
    This script identifies the character from the 1972 film 'Solaris'
    who expresses shame at missing the sound of rustling leaves.
    """

    question = "In the 1972 Andrei Tarkovsky movie Solaris which character is ashamed to miss the sound of leaves rustling on Earth?"

    answer_choices = {
        'A': 'Kris',
        'B': 'Hari',
        'C': 'Snaut',
        'D': 'Sartorius',
        'E': 'Gibarian'
    }

    # Based on the film's dialogue, Dr. Snaut is the character who says this.
    # During a philosophical monologue, he laments humanity's focus on conquering space
    # while neglecting fundamental human feelings and experiences connected to Earth.
    correct_character_name = "Snaut"
    
    # Find the corresponding letter for the correct answer.
    correct_answer_letter = None
    for letter, name in answer_choices.items():
        if name == correct_character_name:
            correct_answer_letter = letter
            break

    # Print the explanation and the result.
    print(f"Movie: Solaris (1972)")
    print(f"Question: {question}")
    print("-" * 20)
    print("Explanation:")
    print("In a key scene, Dr. Snaut delivers a monologue about humanity's purpose in space.")
    print("He expresses a deep sense of longing for Earth and its simple, natural sensations.")
    print("He explicitly states that he is ashamed to be missing the 'rustle of leaves' on Earth.")
    print("-" * 20)
    print(f"The correct character is {correct_character_name}, which corresponds to answer choice {correct_answer_letter}.")

solve_solaris_question()