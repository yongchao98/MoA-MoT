def display_earliest_date():
    """
    This function displays the earliest known date recorded in a pre-Columbian
    aboriginal writing system, which is found on Stela 2 from Chiapa de Corzo, Mexico.
    """
    
    # The date is recorded in the Mesoamerican Long Count system.
    # The date is 7.16.3.2.13.
    baktun = 7
    katun = 16
    tun = 3
    uinal = 2
    kin = 13
    
    print("The earliest known, fully recorded date from a pre-Columbian civilization is found on Stela 2 from Chiapa de Corzo.")
    print("The date is recorded using the Mesoamerican Long Count calendar.")
    print("\n--- Long Count Date Breakdown ---")
    print(f"{baktun:<10} (Baktun cycles of 144,000 days)")
    print(f"{katun:<10} (Katun cycles of 7,200 days)")
    print(f"{tun:<10} (Tun cycles of 360 days)")
    print(f"{uinal:<10} (Uinal cycles of 20 days)")
    print(f"{kin:<10} (Kin, or single days)")
    
    print("\nThe full Long Count date as written is:")
    print(f"{baktun}.{katun}.{tun}.{uinal}.{kin}")
    
    print("\nThis date corresponds to the following in our modern Gregorian calendar:")
    print("December 7, 36 BCE")

display_earliest_date()