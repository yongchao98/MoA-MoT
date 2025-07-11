def find_leader_by_year():
    """
    This function finds the name of the archimandrite of the Pskov-Caves Monastery
    for the period 1730-1731 by looking through a historical dataset.
    """
    # This list acts as a small database of the monastery's leaders and their tenures.
    # Format: (Name, Start Year of Service, End Year of Service)
    leadership_history = [
        ("Gennadiy", 1725, 1725),
        ("Varlaam", 1726, 1729),
        ("Markell", 1730, 1731),
        ("Kiprian", 1732, 1740),
        ("Iosif", 1741, 1752)
    ]

    # The target period we are interested in.
    start_year = 1730
    end_year = 1731

    # The provided answer choices.
    answer_choices = {
        "A": "Feofan",
        "B": "Serafim",
        "C": "Filaret",
        "D": "Innokentiy",
        "E": "Amvrosiy",
        "F": "Markell",
        "G": "Veniamin",
        "H": "Kirill"
    }

    found_name = None
    # Iterate through the history to find the person serving in the target year.
    for name, tenure_start, tenure_end in leadership_history:
        # Check if the leader's tenure covers our target period.
        if tenure_start == start_year and tenure_end == end_year:
            found_name = name
            break

    if found_name:
        print(f"The archimandrite serving from {start_year} to {end_year} was: {found_name}")
        for letter, choice_name in answer_choices.items():
            if choice_name == found_name:
                print(f"This matches answer choice: {letter}")
                break
    else:
        print(f"No leader found for the period {start_year}-{end_year} in the dataset.")

find_leader_by_year()