import textwrap

def find_shared_characters():
    """
    Finds which Shakespearean title characters from a given list
    are also mentioned by name in Dante's 'The Divine Comedy'.
    """
    # Step 1: List all unique Shakespearean title characters from the answer choices.
    shakespeare_title_characters = sorted(list({
        'Julius Caesar', 'Pericles', 'Cleopatra', 'King John', 'Troilus', 'Antony'
    }))

    # Step 2: Create a data source of characters mentioned by name in The Divine Comedy.
    # This is based on literary analysis of the text.
    dante_mentions = {
        'Julius Caesar': 'Inferno, Canto IV: Found in Limbo among the virtuous pagans.',
        'Cleopatra': 'Inferno, Canto V: Found in the Second Circle of Hell with the lustful.',
        # Note: Other characters like Pericles, King John, Troilus, and Antony are not
        # mentioned by name in The Divine Comedy.
    }

    print("Checking for Shakespearean title characters in Dante's The Divine Comedy...")
    print("-" * 70)

    found_characters = []

    # Step 3: Cross-reference the lists.
    for character in shakespeare_title_characters:
        if character in dante_mentions:
            status = "FOUND"
            details = dante_mentions[character]
            found_characters.append(character)
        else:
            status = "NOT FOUND"
            details = f"{character} is not mentioned by name in The Divine Comedy."
        
        print(f"Character: {character:<15} | Status: {status:<10} | Details: {details}")

    print("-" * 70)
    
    # Step 4: Output the final list of matching characters.
    final_result = " and ".join(found_characters)
    print(f"The title characters from the options who are mentioned by name in The Divine Comedy are: {final_result}.")
    print("\nComparing this result with the answer choices, the correct option is the one that lists only Julius Caesar and Cleopatra.")


find_shared_characters()