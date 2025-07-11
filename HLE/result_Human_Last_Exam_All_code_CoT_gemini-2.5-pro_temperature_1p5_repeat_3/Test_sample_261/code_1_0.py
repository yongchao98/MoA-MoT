def find_archimandrite():
    """
    This function finds and prints the name of the archimandrite of the
    Pskov-Caves Monastery for the period 1730-1731.
    """
    # Historical data of the leaders (archimandrites) of the Pskov-Caves Monastery
    # around the specified time. The format is {Name: (start_year, end_year)}.
    abbots_tenure = {
        "Varlaam": (1725, 1729),
        "Markell": (1729, 1731),
        "Innokentiy": (1731, 1740)
    }

    target_start = 1730
    target_end = 1731
    
    found_abbot = None

    # Search for the abbot whose tenure includes the target period
    for abbot, (start_year, end_year) in abbots_tenure.items():
        if start_year <= target_start and end_year >= target_end:
            found_abbot = abbot
            break
            
    if found_abbot:
        print(f"The archimandrite of the Pskov-Caves Monastery from {target_start} to {target_end} was: {found_abbot}")
    else:
        print("Could not find the archimandrite for the specified period in the dataset.")

find_archimandrite()