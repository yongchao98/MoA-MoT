# To run this code, you may first need to install the `lunardate` library by running:
# pip install lunardate

from lunardate import SolarDate

def find_next_matching_birthday():
    """
    Finds the next year after 1980 where the Lunar date for Oct 1st
    matches the Lunar date for 1980-10-01.
    """
    # The person's solar birth date
    birth_year = 1980
    birth_month = 10
    birth_day = 1

    # 1. Convert the solar birth date to its corresponding lunar date
    try:
        solar_birth_date = SolarDate(birth_year, birth_month, birth_day)
        lunar_birth_date = solar_birth_date.to_lunar()
    except ImportError:
        print("Error: The 'lunardate' library is not installed.")
        print("Please install it by running: pip install lunardate")
        return
    except Exception as e:
        print(f"An error occurred: {e}")
        return
        
    target_lunar_month = lunar_birth_date.month
    target_lunar_day = lunar_birth_date.day

    # 2. Start searching from the year after birth
    start_year = birth_year + 1

    # We will search for the next 100 years, which is a safe range
    for year in range(start_year, start_year + 100):
        # 3. Get the lunar date for Oct 1st of the current year
        current_solar_date = SolarDate(year, birth_month, birth_day)
        current_lunar_date = current_solar_date.to_lunar()
        
        # 4. Check if the lunar month and day match the person's lunar birthday
        if (current_lunar_date.month == target_lunar_month and 
            current_lunar_date.day == target_lunar_day):
            # 5. Found the year. Print it and stop.
            print(year)
            return

    print("No matching year found in the next 100 years.")

if __name__ == '__main__':
    find_next_matching_birthday()
