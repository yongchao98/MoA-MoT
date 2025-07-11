def find_character_in_dog_s_heart():
    """
    This script identifies the correct character from a multiple-choice question
    about Mikhail Bulgakov's "A Dog's Heart".
    """
    question = "Who was the older woman Polygraf attempted to assault?"
    
    choices = {
        'A': 'Vasnetsova',
        'B': 'Varvana',
        'C': 'Zina',
        'D': 'Maria',
        'E': 'Darya'
    }
    
    correct_answer_key = 'A'
    correct_answer_name = choices[correct_answer_key]
    
    explanation = (
        "In the novel, Sharikov brings a typist named Vasnetsova to the professor's apartment. "
        "He attempts to coerce her into a relationship and force her to live with him through lies and intimidation. "
        "This incident is the most prominent example of such an assaultive attempt by Sharikov against a woman."
    )
    
    print(f"Question: {question}")
    print("Answer Choices:")
    for key, value in choices.items():
        print(f"{key}. {value}")
    
    print("\n--- SOLUTION ---")
    print(f"The correct choice is {correct_answer_key}: {correct_answer_name}.")
    print("\nExplanation:")
    print(explanation)

find_character_in_dog_s_heart()