import sys

try:
    from lunardate import LunarDate, SolarDate
except ImportError:
    print("Error: The 'lunardate' library is required to run this script.")
    print("Please install it by running: pip install lunardate")
    sys.exit(1)

def find_next_matching_birthday():
    """
    Finds the next year after 1980 where the Lunar birthday corresponding to 
    1980-10-01 falls on the Solar date of October 1st.
    """
    # Define the person's solar birth date
    birth_year = 1980
    birth_month = 10
    birth_day = 1

    # Convert the solar birth date to its lunar equivalent
    try:
        birth_lunar_date = LunarDate.fromSolarDate(birth_year, birth_month, birth_day)
        target_lunar_month = birth_lunar_date.month
        target_lunar_day = birth_lunar_date.day
    except ValueError as e:
        print(f"Error: Could not process the initial birth date. {e}")
        return

    print(f"The person's Solar birthday is on {birth_month}-{birth_day}.")
    print(f"Their corresponding Lunar birthday is Month: {target_lunar_month}, Day: {target_lunar_day}.")
    print("Searching for the next year where their Lunar birthday falls on October 1st...")
    
    # Start searching from the year after birth
    current_year_to_check = birth_year + 1

    while True:
        try:
            # For the current year, create a LunarDate object for their birthday
            lunar_birthday_in_year = LunarDate(current_year_to_check, target_lunar_month, target_lunar_day)
            
            # Convert this lunar date back to its solar equivalent
            solar_equivalent = lunar_birthday_in_year.toSolarDate()

            # Check if the solar equivalent is October 1st
            if solar_equivalent.month == 10 and solar_equivalent.day == 1:
                result_year = solar_equivalent.year
                print(f"\nMatch found!")
                print(f"The next year this occurs is: {result_year}")
                break
                
        except ValueError:
            # This handles cases where the lunar date is invalid for a given year
            # (e.g., day 30 in a 29-day month, or leap month considerations).
            # We simply skip to the next year.
            pass

        current_year_to_check += 1
        
        # A safety break to prevent an unexpected infinite loop
        if current_year_to_check > birth_year + 150:
            print("Could not find a matching year within 150 years.")
            break

if __name__ == "__main__":
    find_next_matching_birthday()