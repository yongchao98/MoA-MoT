def find_bryman_characteristics():
    """
    Identifies the correct pair of Disneyfication characteristics
    discussed by Alan Bryman from a list of choices.
    """
    # Alan Bryman's four main dimensions of Disneyization
    bryman_dimensions = {
        "theming",
        "hybrid consumption",
        "merchandising",
        "performative labor"
    }

    # The answer choices provided to the user
    choices = {
        "A": ["hybrid consumption", "merchandising"],
        "B": ["performative labor", "sanitization"],
        "C": ["trivialization", "theming"],
        "D": ["sanitization", "trivialization"],
        "E": ["disneyfication", "disneyization"],
        "F": ["mcdonaldization", "disneyization"],
        "G": ["theming", "performative labor"]
    }

    print("Analyzing the options based on Alan Bryman's four dimensions of Disneyization...")
    print(f"The core dimensions are: {', '.join(sorted(list(bryman_dimensions)))}.\n")

    correct_options = []
    for letter, characteristics in choices.items():
        char1, char2 = characteristics
        # Check if both items in the choice are part of Bryman's core dimensions
        if char1 in bryman_dimensions and char2 in bryman_dimensions:
            correct_options.append(letter)
            print(f"Option {letter}: '{char1}' and '{char2}'.  -> Found a valid pair.")
        else:
            print(f"Option {letter}: '{char1}' and '{char2}'.  -> Not a valid pair.")

    # Both A and G are technically correct. However, theming is the foundational
    # element, and performative labor is a key sociological concept highlighted by Bryman.
    # Therefore, G is often considered the 'best' answer.
    final_answer = "G"

    print(f"\nMultiple options ({', '.join(correct_options)}) contain two valid characteristics.")
    print(f"However, the most representative pair, including the foundational 'theming' and the sociological 'performative labor', is option {final_answer}.")

find_bryman_characteristics()