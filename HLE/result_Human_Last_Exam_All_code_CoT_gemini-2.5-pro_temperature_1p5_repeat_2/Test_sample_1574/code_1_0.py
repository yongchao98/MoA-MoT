def find_character():
    """
    This function simulates a search through a knowledge base on "A Dog's Heart"
    to find the character who was the target of Sharikov's assault.
    """
    character_interactions = {
        'Zina': "A young maid whom Sharikov harasses and makes inappropriate advances towards within the apartment.",
        'Darya': "The cook, towards whom Sharikov is rude and threatening.",
        'Vasnetsova': "A typist at Sharikov's workplace. Sharikov attempts to coerce her into a relationship by threatening to have her fired. This abuse of power is a key plot point."
    }

    # The question concerns an attempted assault, which strongly points to
    # the incident of coercion and abuse of power, a central theme of the story.
    target_description_keyword = "coerce"

    for name, description in character_interactions.items():
        if target_description_keyword in description:
            print(f"The woman Polygraf Polygrafovich Sharikov attempted to assault by coercing her and threatening her job was: {name}")
            return

find_character()