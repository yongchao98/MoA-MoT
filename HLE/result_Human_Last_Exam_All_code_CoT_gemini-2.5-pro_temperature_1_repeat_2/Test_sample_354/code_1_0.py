def find_bryman_characteristics():
    """
    Identifies the correct pair of Disneyfication characteristics
    from a list of options, based on Alan Bryman's work.
    """
    # Alan Bryman's four key dimensions of Disneyization
    bryman_dimensions = {
        "theming",
        "hybrid consumption",
        "merchandising",
        "performative labor"
    }

    # The multiple-choice options provided
    options = {
        "A": {"hybrid consumption", "merchandising"},
        "B": {"performative labor", "sanitization"},
        "C": {"trivialization", "theming"},
        "D": {"sanitization", "trivialization"},
        "E": {"disneyfication", "disneyization"},
        "F": {"mcdonaldization", "disneyization"},
        "G": {"theming", "performative labor"}
    }

    correct_answer_key = None
    correct_answer_values = None

    # Find the option where both items are in Bryman's list of dimensions
    for key, values in options.items():
        if values.issubset(bryman_dimensions):
            correct_answer_key = key
            correct_answer_values = values
            # There might be more than one correct option (A and G),
            # but we will return the first one found.
            break

    print("In 'The Disneyization of Society', Alan Bryman identifies four main characteristics:")
    for dim in sorted(list(bryman_dimensions)):
        print(f"- {dim.title()}")

    print("\nChecking the options...")
    if correct_answer_key:
        value_list = sorted(list(correct_answer_values))
        print(f"Option {correct_answer_key} ('{value_list[0]}' and '{value_list[1]}') correctly lists two of these characteristics.")
        print("\n<<<" + correct_answer_key + ">>>")
    else:
        print("Could not programmatically determine the correct answer from the provided options.")

# Execute the function to find and print the answer
find_bryman_characteristics()