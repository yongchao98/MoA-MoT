def solve():
    """
    Calculates and displays the earliest recorded date from a pre-Columbian
    civilization, found on Stela 2 at Chiapa de Corzo.
    """
    print("The earliest known recorded date is from Stela 2, Chiapa de Corzo, Mexico.")
    print("The date is recorded in the Maya Long Count system as: 7.16.3.2.13\n")

    # 1. Define the Long Count date components
    baktun = 7
    katun = 16
    tun = 3
    uinal = 2
    kin = 13

    # Days in each Long Count period
    days_in_baktun = 144000
    days_in_katun = 7200
    days_in_tun = 360
    days_in_uinal = 20
    days_in_kin = 1

    # 2. Calculate the total number of days from the epoch start
    print("Calculating the total number of days from the calendar's start point:")
    
    baktun_days = baktun * days_in_baktun
    print(f"{baktun} Baktun * {days_in_baktun} days = {baktun_days} days")

    katun_days = katun * days_in_katun
    print(f"{katun} Katun * {days_in_katun} days = {katun_days} days")

    tun_days = tun * days_in_tun
    print(f"{tun} Tun * {days_in_tun} days = {tun_days} days")
    
    uinal_days = uinal * days_in_uinal
    print(f"{uinal} Uinal * {days_in_uinal} days = {uinal_days} days")
    
    kin_days = kin * days_in_kin
    print(f"{kin} Kin * {days_in_kin} day = {kin_days} days")

    total_days = baktun_days + katun_days + tun_days + uinal_days + kin_days
    print(f"\nTotal days since epoch: {baktun_days} + {katun_days} + {tun_days} + {uinal_days} + {kin_days} = {total_days} days\n")

    # 3. Convert to Gregorian date using Julian Day Number (JDN)
    # The JDN for the Maya epoch start (0.0.0.0.0) using the GMT correlation is 584285.
    mayan_epoch_jdn = 584285
    print(f"The Julian Day Number (JDN) for the Maya epoch start is: {mayan_epoch_jdn}")

    target_jdn = mayan_epoch_jdn + total_days
    print(f"The JDN for the target date is: {mayan_epoch_jdn} + {total_days} = {target_jdn}\n")
    
    # 4. Function to convert JDN to Gregorian Date
    def jd_to_gregorian(jd):
        j = jd + 32044
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
        return year, month, day

    year, month, day = jd_to_gregorian(target_jdn)
    
    # Format the final date string
    month_names = ["January", "February", "March", "April", "May", "June", 
                   "July", "August", "September", "October", "November", "December"]
    
    # Convert astronomical year (which includes year 0) to BCE
    if year < 1:
        year_str = f"{-year + 1} BCE"
    else:
        year_str = f"{year} CE"

    final_date = f"{month_names[month-1]} {day}, {year_str}"
    
    print("This corresponds to the Gregorian calendar date:")
    print(final_date)
    
    # The prompt requests the final answer in a specific format at the end.
    # The question is "What is the earliest known date...?"
    # The answer is the date itself.
    global final_answer_value
    final_answer_value = final_date

# Run the solver
solve()
<<<December 8, 36 BCE>>>