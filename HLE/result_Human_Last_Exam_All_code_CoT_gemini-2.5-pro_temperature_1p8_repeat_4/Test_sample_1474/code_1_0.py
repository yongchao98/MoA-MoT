import sys

def find_common_theme():
    """
    Analyzes the works of directors Fritz Lang and William Friedkin
    to find a common theme among the given choices.
    """
    # Knowledge base of themes in the directors' films.
    # This is based on film history analysis.
    knowledge_base = {
        "Fritz Lang": {
            "Aboriginal masks": False,
            "Magic wands": False,
            "The first ever cyborgs on screen": True, # For "Metropolis"
            "Bugs": True # For "Die Spinnen" (The Spiders) and themes of infestation/paranoia
        },
        "William Friedkin": {
            "Aboriginal masks": False,
            "Magic wands": False,
            "The first ever cyborgs on screen": False,
            "Bugs": True # For "The Exorcist" (Pazuzu/locusts) and the film "Bug"
        }
    }

    choices = {
        "A": "Aboriginal masks",
        "B": "Magic wands",
        "C": "The first ever cyborgs on screen",
        "D": "Bugs"
    }

    correct_answer_key = "E"
    final_reasoning = "No common theme found from choices A-D."

    print("Analyzing themes for directors Fritz Lang and William Friedkin...\n")

    # Iterate through the choices to find a commonality
    for key, theme in choices.items():
        is_in_lang = knowledge_base["Fritz Lang"][theme]
        is_in_friedkin = knowledge_base["William Friedkin"][theme]

        print(f"Checking for theme: '{theme}'")
        print(f"-> Present in Lang's work: {is_in_lang}")
        print(f"-> Present in Friedkin's work: {is_in_friedkin}")

        if is_in_lang and is_in_friedkin:
            correct_answer_key = key
            # This serves as the "equation" showing the components of the final answer
            final_reasoning = f"The common theme is '{theme}'.\nFritz Lang's work ('The Spiders') and William Friedkin's work ('The Exorcist', 'Bug') both utilize this theme."
            print("Conclusion: Found a common theme!\n")
            break
        else:
            print("Conclusion: Not a common theme.\n")

    print("--- FINAL RESULT ---")
    print(final_reasoning)
    print(f"The correct option is {correct_answer_key}.")
    
    # Required final answer format for the platform
    # Note: The platform requires this specific format. 
    # In a real script, this would be handled differently.
    sys.stdout.write(f"\n<<<{correct_answer_key}>>>")


find_common_theme()