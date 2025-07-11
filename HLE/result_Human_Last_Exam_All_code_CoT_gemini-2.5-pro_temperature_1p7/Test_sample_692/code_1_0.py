def find_earliest_date():
    """
    This function provides information about the earliest known date
    recorded by a pre-Columbian civilization in the Americas.
    """
    
    # The components of the Mesoamerican Long Count date from Chiapa de Corzo Stela 2
    baktun = 7
    katun = 16
    tun = 3
    uinal = 2
    kin = 13

    # The corresponding date in the proleptic Gregorian calendar
    gregorian_date = "December 7, 36 BC"
    
    # The artifact and location
    artifact = "Stela 2"
    location = "Chiapa de Corzo, Mexico"
    
    print("The earliest known recorded date from a pre-Columbian civilization is found on {} from {}.".format(artifact, location))
    print("In the Long Count calendar, the date is written as:")
    print(f"{baktun}.{katun}.{tun}.{uinal}.{kin}")
    print(f"This corresponds to the date: {gregorian_date}.")

if __name__ == '__main__':
    find_earliest_date()