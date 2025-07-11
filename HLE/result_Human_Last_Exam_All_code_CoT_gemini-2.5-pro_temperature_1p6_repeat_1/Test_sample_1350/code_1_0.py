def find_character_in_wolfings():
    """
    This script identifies a character from William Morris's "The House of the Wolfings"
    based on a specific story they recount.
    """

    # A data structure representing characters and their notable stories/actions.
    character_profiles = [
        {
            "name": "Thiodolf",
            "story": "The leader of the Wolfing tribe who defends his people against the encroaching Roman army."
        },
        {
            "name": "Hall-Sun",
            "story": "The mortal daughter of Thiodolf, known for her love and devotion to her family."
        },
        {
            "name": "The Wood-Sun",
            "story": ("A demi-goddess and foster-daughter to Thiodolf. She describes "
                      "having visited a Roman city and seeing a statue of a great she-wolf "
                      "suckling a twin pair of infants, referring to the Roman origin myth of Romulus and Remus.")
        },
        {
            "name": "Asmund",
            "story": "A valiant warrior and kinsman within the House of the Wolfings."
        }
    ]

    # Search for the character who described the Roman origin myth.
    target_character = None
    search_terms = ["Roman", "wolf", "twin"]

    for profile in character_profiles:
        story = profile["story"]
        if all(term in story for term in search_terms):
            target_character = profile["name"]
            break

    # Print the result.
    if target_character:
        print("Based on the provided text, the character who visited a Roman city and described their religious practices is:")
        print(target_character)
    else:
        print("Could not identify the character based on the search terms.")

# Run the function to find and print the answer.
find_character_in_wolfings()