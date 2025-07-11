import collections

def find_common_theme():
    """
    Analyzes the filmographies of Fritz Lang and William Friedkin
    to find a common theme from a given list of choices.
    """
    # Step 1: Create a knowledge base of themes for each director.
    # This includes literal themes and common critical interpretations.
    knowledge_base = {
        "Fritz Lang": {
            "The first ever cyborgs on screen (in 'Metropolis')",
            "Metaphorical bugs/insects (the swarming, dehumanized workers in 'Metropolis')",
            "Paranoia and surveillance (in 'Dr. Mabuse the Gambler' and 'M')"
        },
        "William Friedkin": {
            "Bugs/insects (the demon Pazuzu is linked to locusts in 'The Exorcist'; the central theme of 'Bug')",
            "Unexplained mysteries and faith crisis (in 'The Exorcist')",
            "Moral ambiguity and crime (in 'The French Connection')"
        }
    }

    # Step 2: Define the answer choices and their corresponding keywords.
    answer_choices = {
        "A": "Aboriginal masks",
        "B": "Magic wands",
        "C": "cyborgs",
        "D": "Bugs"
    }

    print("Analyzing cinematic themes for Fritz Lang and William Friedkin...\n")

    correct_answer_letter = None
    correct_answer_text = ""

    # Step 3: Iterate through choices and check for common themes.
    for letter, choice_text in answer_choices.items():
        keyword = choice_text.lower().split(" ")[-1] # a simple way to get the main keyword, e.g., 'bugs'

        lang_has_theme = False
        for theme in knowledge_base["Fritz Lang"]:
            if keyword in theme.lower():
                lang_has_theme = True
                break
        
        friedkin_has_theme = False
        for theme in knowledge_base["William Friedkin"]:
            if keyword in theme.lower():
                friedkin_has_theme = True
                break

        # Check if the theme is present for both directors
        is_common_theme = lang_has_theme and friedkin_has_theme
        
        print(f"Checking Option {letter}: '{choice_text}'")
        print(f"  - Theme found in Lang's work: {lang_has_theme}")
        print(f"  - Theme found in Friedkin's work: {friedkin_has_theme}")
        print(f"  - Is it a common theme? {'Yes' if is_common_theme else 'No'}\n")

        if is_common_theme:
            correct_answer_letter = letter
            correct_answer_text = choice_text

    # Step 4: Print the final conclusion.
    if correct_answer_letter:
        print("Conclusion: The correct answer is the theme found in both oeuvres.")
        print(f"Final Answer: {correct_answer_letter}. {correct_answer_text}")
    else:
        print("Conclusion: No single matching theme found for both directors based on the provided choices.")


if __name__ == "__main__":
    find_common_theme()