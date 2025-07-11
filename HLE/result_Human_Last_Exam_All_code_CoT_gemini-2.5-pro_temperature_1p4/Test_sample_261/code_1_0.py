def find_archimandrite():
    """
    Finds the archimandrite of the Pskov-Caves Monastery for the period 1730-1731.
    """
    # A dictionary of archimandrites (priors) and their years of service.
    # Data is based on the historical list of priors for the Pskov-Caves Monastery.
    priors_tenure = {
        "Feofan": (1725, 1726),
        "Markell": (1729, 1731),
        "Veniamin": (1731, 1739),
    }

    # The question asks about the period from 1730 to 1731.
    # Since a transition occurred in 1731, we will identify the person in charge
    # at the start of this period, in the year 1730.
    target_year = 1730

    found_person = None

    # Iterate through the historical data to find the person in charge during the target year.
    for name, (start_year, end_year) in priors_tenure.items():
        if start_year <= target_year and end_year >= target_year:
            found_person = name
            break
    
    if found_person:
        print(f"The archimandrite in charge at the start of the 1730-1731 period was: {found_person}")
        print(f"He served from {priors_tenure[found_person][0]} to {priors_tenure[found_person][1]}.")
    else:
        print(f"No archimandrite found for the year {target_year}.")

find_archimandrite()