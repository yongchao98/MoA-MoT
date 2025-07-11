def find_common_element():
    """
    Analyzes cinematic elements in the works of Fritz Lang and William Friedkin
    to find a common theme from a given list of choices.
    """
    
    # A knowledge base representing key elements in each director's filmography.
    director_facts = {
        "Fritz Lang": {
            "Cyborgs": "Yes, 'Metropolis' (1927) features the Maschinenmensch, considered the first major cyborg/robot in cinema.",
            "Bugs": "Yes, a character in 'Dr. Mabuse the Gambler' (1922) suffers a hallucination of swarming insects, a potent image of psychological collapse."
        },
        "William Friedkin": {
            "Pazuzu Statue": "Yes, 'The Exorcist' (1973) opens with the discovery of a Pazuzu artifact, which is Mesopotamian, not Aboriginal.",
            "Bugs": "Yes, the film 'Bug' (2006) is centered on paranoia about insect infestation, and 'The Exorcist' features the iconic, insect-like 'spider-walk' scene."
        }
    }

    # The multiple-choice options provided in the question.
    answer_choices = {
        "A": "Aboriginal masks",
        "B": "Magic wands",
        "C": "The first ever cyborgs on screen",
        "D": "Bugs",
        "E": "None of the above"
    }
    
    # We will map our knowledge base keys to the options.
    keyword_map = {
        "A": None, # 'Aboriginal masks' is distinct from 'Pazuzu Statue'.
        "B": None, # No specific data on 'Magic wands'.
        "C": "Cyborgs",
        "D": "Bugs"
    }

    correct_choice = "E"  # Default assumption
    
    print("Analyzing the filmographies of Fritz Lang and William Friedkin...\n")

    for choice_letter, choice_text in answer_choices.items():
        if choice_letter == "E":
            continue

        print(f"Checking Option {choice_letter}: '{choice_text}'")
        
        keyword = keyword_map.get(choice_letter)
        if not keyword:
            print("-> This element is not a noted feature in the works of either director.\n")
            continue
            
        is_in_lang = keyword in director_facts["Fritz Lang"]
        is_in_friedkin = keyword in director_facts["William Friedkin"]

        if is_in_lang:
            print(f" - Lang connection: {director_facts['Fritz Lang'][keyword]}")
        else:
            print(" - Lang connection: Not found.")

        if is_in_friedkin:
            print(f" - Friedkin connection: {director_facts['William Friedkin'][keyword]}")
        else:
            print(" - Friedkin connection: Not found.")
            
        if is_in_lang and is_in_friedkin:
            correct_choice = choice_letter
            print(f"-> Conclusion: '{choice_text}' is a common element.\n")
            break  # Found the answer
        else:
            print(f"-> Conclusion: '{choice_text}' is not found in the work of both directors.\n")
            
    print("--- FINAL RESULT ---")
    print(f"The common element shared by both Fritz Lang and William Friedkin is: '{answer_choices[correct_choice]}'.")
    print(f"The correct option is: {correct_choice}")

find_common_element()