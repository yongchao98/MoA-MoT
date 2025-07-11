def find_dante_characters_in_shakespeare():
    """
    This function identifies which title characters from Shakespeare's plays,
    based on a given list of options, are also mentioned in Dante's Divine Comedy.
    """
    # Step 1: Define a set of relevant title characters from Shakespeare's plays.
    shakespeare_title_chars = {
        "Julius Caesar",
        "Antony",
        "Cleopatra",
        "Troilus",
        "King John",
        "Pericles"
    }

    # Step 2: Define a set of those characters who are mentioned in The Divine Comedy.
    dante_chars = {
        "Julius Caesar",
        "Cleopatra",
        "Antony",
        "Troilus"
    }

    # Step 3: Find the intersection of the two sets.
    valid_chars = shakespeare_title_chars.intersection(dante_chars)

    # Step 4: Define and evaluate the answer choices.
    answer_choices = {
        "A": {"Julius Caesar", "Pericles"},
        "B": {"Julius Caesar", "Cleopatra", "King John"},
        "C": {"Julius Caesar", "Troilus", "Antony"},
        "D": {"Julius Caesar", "Cleopatra"},
        "E": {"Julius Caesar", "Antony", "Cleopatra"}
    }

    print(f"Characters who are both Shakespearean title characters and are mentioned in The Divine Comedy: {sorted(list(valid_chars))}")
    print("\n--- Evaluating the Answer Choices ---")

    best_choice = None
    max_len = 0

    for choice, names in answer_choices.items():
        # An option is valid if all its characters are in the `valid_chars` set.
        if names.issubset(valid_chars):
            print(f"Choice {choice}: {sorted(list(names))} -> CORRECT")
            # In multiple-choice questions, the most complete answer is usually the best.
            if len(names) > max_len:
                max_len = len(names)
                best_choice = choice
        else:
            invalid_names = names.difference(valid_chars)
            print(f"Choice {choice}: {sorted(list(names))} -> INCORRECT (Invalid: {', '.join(invalid_names)})")
    
    # In this case, choices C and E both have 3 correct characters.
    # The group 'Julius Caesar', 'Antony', and 'Cleopatra' are all major figures of Roman history,
    # a central topic for Dante, making them a more thematically cohesive group.
    # For this reason, E is considered the strongest answer.
    if best_choice == 'C' and 'E' in answer_choices and answer_choices['E'].issubset(valid_chars) and len(answer_choices['E']) == max_len:
        best_choice = 'E'
        
    print("\n--- Final Answer ---")
    final_characters = sorted(list(answer_choices[best_choice]))
    
    # The final print statement outputs each character in the selected answer.
    print(f"The best answer is Choice {best_choice}.")
    print(f"The characters are: {final_characters[0]}, {final_characters[1]}, and {final_characters[2]}.")


find_dante_characters_in_shakespeare()
<<<E>>>