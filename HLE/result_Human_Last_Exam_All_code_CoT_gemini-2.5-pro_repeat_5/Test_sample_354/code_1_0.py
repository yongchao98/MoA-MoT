def solve_bryman_question():
    """
    Analyzes the multiple-choice question about Alan Bryman's
    characteristics of Disneyfication.
    """
    # Step 1: Define the four core characteristics identified by Bryman.
    bryman_characteristics = {
        "theming",
        "hybrid consumption",
        "merchandising",
        "performative labor"
    }

    # Step 2: Define the answer choices provided.
    options = {
        'A': {"hybrid consumption", "merchandising"},
        'B': {"performative labor", "sanitization"},
        'C': {"trivialization", "theming"},
        'D': {"sanitization", "trivialization"},
        'E': {"Disneyfication", "Disneyization"},
        'F': {"McDonaldization", "Disneyization"},
        'G': {"theming", "performative labor"}
    }

    print("Evaluating the options against Bryman's four dimensions of Disneyization:")
    
    # Step 3 & 4: Check each option and identify the valid ones.
    valid_options = []
    for key, value_set in options.items():
        # A choice is valid if its characteristics are a subset of Bryman's.
        if value_set.issubset(bryman_characteristics):
            valid_options.append(key)
            print(f"- Option {key}: Contains a valid pair of characteristics: {', '.join(value_set)}.")
        else:
            # Find which parts are not in Bryman's core list for clarity.
            invalid_parts = value_set.difference(bryman_characteristics)
            if invalid_parts:
                print(f"- Option {key}: Contains an invalid characteristic: {', '.join(invalid_parts)}.")
            else: # This case is for options like E and F which are names of concepts.
                 print(f"- Option {key}: Does not represent two characteristics of the process.")


    # Step 5: Explain the final choice if ambiguity exists.
    print("\n--- Conclusion ---")
    if len(valid_options) > 1:
        print(f"Analysis shows that options {', '.join(valid_options)} both contain correct pairs of characteristics.")
        print("However, 'theming' (creating a narrative environment) and 'performative labor' (employees acting a role) are often considered the primary pillars that construct the unique, immersive experience of Disneyization.")
        print("Therefore, Option G is the strongest choice.")
        final_answer = "G"
    elif len(valid_options) == 1:
        final_answer = valid_options[0]
        print(f"The only correct option is {final_answer}.")
    else:
        final_answer = "None"
        print("No option correctly lists two of Bryman's characteristics.")

# Execute the analysis function.
solve_bryman_question()