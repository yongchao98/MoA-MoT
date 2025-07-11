def solve_shakespeare_dante_question():
    """
    This script determines which Shakespearean title characters from a given list
    are also mentioned by name in Dante's "The Divine Comedy".
    """

    # Database of Shakespearean title characters mentioned by name in The Divine Comedy.
    # This is based on literary analysis of Dante's text.
    in_divine_comedy = {"Julius Caesar", "Cleopatra"}

    # The provided answer choices as a dictionary.
    answer_choices = {
        "A": ["Julius Caesar", "Pericles"],
        "B": ["Julius Caesar", "Cleopatra", "King John"],
        "C": ["Julius Caesar", "Troilus", "Antony"],
        "D": ["Julius Caesar", "Cleopatra"],
        "E": ["Julius Caesar", "Antony", "Cleopatra"],
    }

    print("Analyzing which Shakespearean characters appear in Dante's Divine Comedy...")
    print("-" * 60)
    print("Factual Check:")
    print("- Julius Caesar: YES, mentioned in Inferno, Canto IV.")
    print("- Cleopatra: YES, mentioned in Inferno, Canto V.")
    print("- Antony: NO, not mentioned by name.")
    print("- Troilus: NO, not mentioned by name.")
    print("- Pericles: NO, not mentioned.")
    print("- King John: NO, not mentioned.")
    print("-" * 60)

    correct_choice_letter = None
    # Find the choice where all characters listed are in our "database".
    for choice, characters in answer_choices.items():
        if all(character in in_divine_comedy for character in characters):
            # A potential issue could be if one choice is a subset of another correct choice.
            # We must verify if the list is exhaustive. In this case, only one choice is fully correct.
            if len(characters) == len(in_divine_comedy):
                correct_choice_letter = choice
                break
    
    # Fallback in case the choices are subsets
    if not correct_choice_letter:
        for choice, characters in answer_choices.items():
            if all(character in in_divine_comedy for character in characters):
                correct_choice_letter = choice
                break


    print(f"The correct option is {correct_choice_letter}.")
    print("The characters from this option who appear by name in the poem are:")
    final_characters = answer_choices[correct_choice_letter]
    
    # Per instruction, output each name in the final correct answer
    for char in final_characters:
        print(f"-> {char}")

solve_shakespeare_dante_question()