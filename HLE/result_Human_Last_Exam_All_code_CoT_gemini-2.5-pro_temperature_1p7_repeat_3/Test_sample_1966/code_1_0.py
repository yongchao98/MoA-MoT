def solve_shakespeare_dante_question():
    """
    This script identifies which title characters from Shakespeare's plays,
    among the given options, are mentioned by name in Dante's 'The Divine Comedy'.
    """

    # Step 1: Create a set of Shakespearean title characters that are verifiably
    # mentioned by name in 'The Divine Comedy'.
    # - Julius Caesar: Yes, in Inferno, Canto IV.
    # - Cleopatra: Yes, in Inferno, Canto V.
    # - Antony, Pericles, King John, Troilus: No, these characters are not mentioned by name.
    characters_in_divine_comedy = {'Julius Caesar', 'Cleopatra'}

    # Step 2: Define the character lists for each multiple-choice option.
    answer_choices = {
        'A': ['Julius Caesar', 'Pericles'],
        'B': ['Julius Caesar', 'Cleopatra', 'King John'],
        'C': ['Julius Caesar', 'Troilus', 'Antony'],
        'D': ['Julius Caesar', 'Cleopatra'],
        'E': ['Julius Caesar', 'Antony', 'Cleopatra']
    }

    print("Checking which characters from each option are in The Divine Comedy...\n")

    correct_choice = None
    
    # Step 3: Iterate through each option and check for validity.
    for choice, characters in answer_choices.items():
        # Check if all characters in the current choice are present in our confirmed set.
        is_correct_subset = set(characters).issubset(characters_in_divine_comedy)
        
        print(f"--- Option {choice} ---")
        print(f"Characters: {', '.join(characters)}")
        
        # Verify each character individually for a detailed explanation.
        for character in characters:
            if character in characters_in_divine_comedy:
                print(f"-> '{character}' is mentioned in The Divine Comedy.")
            else:
                print(f"-> '{character}' is NOT mentioned in The Divine Comedy.")
        
        print(f"Result for Option {choice}: All characters mentioned? {is_correct_subset}\n")

        if is_correct_subset:
            correct_choice = choice
            
    print("="*30)
    print(f"The correct choice is {correct_choice}.")
    print("The characters are: " + ", ".join(answer_choices[correct_choice]))
    print("="*30)

# Run the solver function
solve_shakespeare_dante_question()