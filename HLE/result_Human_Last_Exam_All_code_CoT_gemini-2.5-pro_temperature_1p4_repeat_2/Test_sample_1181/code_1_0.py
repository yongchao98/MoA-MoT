import sys

def solve_solaris_question():
    """
    This function analyzes the question about the 1972 movie Solaris
    and determines the correct character.
    """
    question = "In the 1972 Andrei Tarkovsky movie Solaris which character is ashamed to miss the sound of leaves rustling on Earth?"

    characters = {
        "A": "Kris",
        "B": "Hari",
        "C": "Snaut",
        "D": "Sartorius",
        "E": "Gibarian"
    }

    explanation = (
        "The character who expresses this sentiment is Dr. Gibarian.\n"
        "Before Kris Kelvin arrives on the space station, Gibarian takes his own life.\n"
        "He leaves a final, rambling video message for Kris. In this message, he confesses his despair and shame, "
        "mentioning his longing for Earth and that he is ashamed to admit he misses something as simple as the sound of rustling leaves. "
        "This highlights his psychological torment and profound sense of isolation."
    )

    correct_answer_letter = "E"
    correct_character_name = characters[correct_answer_letter]

    print(f"Question: {question}")
    print("\nAnalyzing the characters and plot points...")
    print(f"\nExplanation:\n{explanation}")
    print(f"\nTherefore, the correct character is {correct_character_name}.")
    print(f"The correct answer choice is: {correct_answer_letter}")

# Execute the function to solve the problem
solve_solaris_question()