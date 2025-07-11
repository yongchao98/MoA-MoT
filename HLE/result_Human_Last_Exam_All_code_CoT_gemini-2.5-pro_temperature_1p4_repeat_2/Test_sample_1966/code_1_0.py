def find_shared_characters():
    """
    Identifies which title characters from Shakespeare's plays are also mentioned by name in Dante's "The Divine Comedy".
    """
    # Step 1: List all unique Shakespearean title characters from the answer choices.
    shakespearean_characters = [
        "Julius Caesar",
        "Pericles",
        "Cleopatra",
        "King John",
        "Troilus",
        "Antony"
    ]

    # Step 2: Create a dataset indicating which characters are mentioned by name in The Divine Comedy.
    # This information is based on well-established literary analysis of Dante's work.
    # True means the character is explicitly named.
    dante_mentions = {
        "Julius Caesar": {"mentioned": True, "location": "Inferno, Canto IV (Limbo)"},
        "Cleopatra": {"mentioned": True, "location": "Inferno, Canto V (Circle of the Lustful)"},
        "Antony": {"mentioned": False, "location": "Not explicitly named, though associated with Cleopatra."},
        "Troilus": {"mentioned": False, "location": "Not explicitly named."},
        "Pericles": {"mentioned": False, "location": "Not mentioned."},
        "King John": {"mentioned": False, "location": "Not mentioned."}
    }

    # Step 3 & 4: Find the intersection of the two sets.
    shared_characters = []
    print("Checking which Shakespearean title characters are named in Dante's Divine Comedy:\n")
    for char in shakespearean_characters:
        if dante_mentions.get(char) and dante_mentions[char]["mentioned"]:
            shared_characters.append(char)
            print(f"- {char}: Yes, mentioned in {dante_mentions[char]['location']}.")
        else:
            print(f"- {char}: No, not explicitly named in The Divine Comedy.")

    # Step 5 & 6: Print the final result.
    print("\n---------------------------------------------------------")
    print("The title characters from Shakespeare's plays who are also mentioned by name in The Divine Comedy are:")
    # Using ' and ' to join the last two elements for better readability
    if len(shared_characters) > 1:
        result_string = ", ".join(shared_characters[:-1]) + " and " + shared_characters[-1]
    elif shared_characters:
        result_string = shared_characters[0]
    else:
        result_string = "None"
        
    print(result_string)
    print("---------------------------------------------------------")
    print("This corresponds to answer choice D.")

find_shared_characters()