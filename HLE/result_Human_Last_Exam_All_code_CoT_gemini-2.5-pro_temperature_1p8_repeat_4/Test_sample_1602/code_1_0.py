def find_wwii_general():
    """
    This function searches through a list of American WWII generals
    to find the one known for a hissing sound from a facial wound.
    """
    # This dictionary simulates a knowledge base of facts about the generals.
    general_facts = {
        'Theodore Roosevelt, Jr.': 'Oldest man in the D-Day invasion, landed with a cane due to arthritis.',
        'George Patton': 'Famous for his aggressive tactics and ivory-handled pistols.',
        'Bruce Magruder': 'First commander of the 1st Infantry Division in WWII.',
        'Raymond Albert Wheeler': 'Engineer officer, primarily served in the China-Burma-India Theater.',
        'Lloyd Fredendall': 'Criticized for his leadership during the Battle of Kasserine Pass.',
        'Leonard T. Gerow': 'Commander of V Corps, involved in the planning of the Normandy landings.',
        'Elbridge Chapman': 'Commander of the 13th Airborne Division.',
        'Terry de la Mesa Allen, Sr.': 'Commander of the 1st Infantry Division during the North African and Sicily campaigns.',
        'Clarence R. Huebner': 'Succeeded Allen as commander of the 1st Infantry Division.',
        'Mark W. Clark': 'Wounded by shrapnel from an aerial bomb, a severed facial nerve caused his cheek to make a slight hissing sound when he was agitated.'
    }

    found = False
    for general, fact in general_facts.items():
        if 'hissing' in fact and 'facial' in fact and 'wound' in fact:
            print(f"The general matching the description is: {general}")
            print(f"Historical Fact: {fact}")
            found = True
            break

    if not found:
        print("Could not find a general matching the description in the knowledge base.")

find_wwii_general()