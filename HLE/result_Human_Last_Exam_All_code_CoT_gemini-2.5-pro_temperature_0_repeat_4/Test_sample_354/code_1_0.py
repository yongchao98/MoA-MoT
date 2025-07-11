def find_bryman_characteristics():
    """
    This script identifies two characteristics of Disneyfication as discussed by Alan Bryman
    by stating his core theory and matching it to the provided options.
    """
    # In "The Disneyization of Society" (2004), Alan Bryman outlines four key dimensions.
    bryman_dimensions = [
        "theming",
        "hybrid consumption",
        "merchandising",
        "performative labor"
    ]

    print("Alan Bryman's four main dimensions of Disneyization are:")
    for dimension in bryman_dimensions:
        print(f"- {dimension.title()}")

    print("\nWe must find the answer choice that contains two of these four dimensions.")

    # Let's analyze the correct choice.
    # Choice G is "theming and performative labor".
    correct_choice_letter = "G"
    characteristic1 = "theming"
    characteristic2 = "performative labor"

    # Check if the characteristics from choice G are in Bryman's list.
    is_correct = characteristic1 in bryman_dimensions and characteristic2 in bryman_dimensions

    if is_correct:
        print(f"\nChoice {correct_choice_letter} lists '{characteristic1}' and '{characteristic2}'.")
        print("Both of these are core dimensions of Disneyization identified by Bryman.")
        print("Therefore, this is the correct answer.")
    else:
        # This part of the code would run if the hardcoded answer was wrong.
        print("\nCould not verify the hardcoded correct answer.")

find_bryman_characteristics()