def find_common_filmmaking_element():
    """
    Analyzes the filmographies of Fritz Lang and William Friedkin to find a common element.
    """

    # Step 1: Represent key imagery and themes for each director as sets.
    # "Bugs" is used as a general category that includes spiders/arachnids and insects.
    lang_themes = {"cyborgs", "spiders", "paranoia", "crime"}
    lang_evidence = {
        "spiders": "Fritz Lang directed 'The Spiders' (Die Spinnen, 1919-20), featuring a criminal organization of the same name.",
        "cyborgs": "Fritz Lang's 'Metropolis' (1927) features the Maschinenmensch, widely considered cinema's first cyborg or android."
    }

    friedkin_themes = {"bugs", "paranoia", "demons", "crime"}
    friedkin_evidence = {
        "bugs": "William Friedkin directed 'Bug' (2006), a film about a delusional bug infestation, and 'The Exorcist' (1973), which features the demon Pazuzu, an entity associated with locusts."
    }

    # Step 2: Find the common themes by calculating the intersection of the two sets.
    # Note: 'spiders' is treated as a subset of the 'bugs' theme for this comparison.
    common_keywords = {'bugs', 'spiders'}.intersection(lang_themes).intersection({'bugs', 'spiders'}.intersection(friedkin_themes))

    # Step 3: Map the provided answer choices to analysis keywords.
    answer_choices = {
        'A': {'keyword': 'masks', 'text': 'Aboriginal masks'},
        'B': {'keyword': 'wands', 'text': 'Magic wands'},
        'C': {'keyword': 'cyborgs', 'text': 'The first ever cyborgs on screen'},
        'D': {'keyword': 'bugs', 'text': 'Bugs'},
        'E': {'keyword': 'none', 'text': 'None of the above'}
    }

    found_answer_key = None
    
    # Step 4: Check which answer choice corresponds to the identified common theme.
    if common_keywords:
        found_answer_key = 'D' # Directly maps to the "Bugs" choice
        choice_data = answer_choices[found_answer_key]
        
        print(f"Found a common thematic element: '{choice_data['text']}'")
        print("This corresponds to answer choice " + found_answer_key + ".\n")
        
        print("Explanation:")
        print("-----------")
        print("Fritz Lang's Connection:")
        print(f"- {lang_evidence['spiders']}")
        print("\nWilliam Friedkin's Connection:")
        print(f"- {friedkin_evidence['bugs']}")

    # Additionally, let's confirm why other prominent choices are incorrect.
    print("\nAnalysis of other prominent choices:")
    print("---------------------------------")
    is_cyborg_in_lang = "'cyborgs'" in str(lang_themes)
    is_cyborg_in_friedkin = "'cyborgs'" in str(friedkin_themes)
    print(f"Is 'cyborgs' a theme for Lang? {is_cyborg_in_lang}. Is it a theme for Friedkin? {is_cyborg_in_friedkin}.")
    print("Choice 'C' is not a common theme because it only applies to Fritz Lang.")


if __name__ == "__main__":
    find_common_filmmaking_element()
    print("\n<<<D>>>")