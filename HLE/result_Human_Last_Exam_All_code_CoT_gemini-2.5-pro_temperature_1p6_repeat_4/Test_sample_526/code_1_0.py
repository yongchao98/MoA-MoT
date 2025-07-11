# First, you might need to install the required library:
# pip install lunardate

from lunardate import SolarDate

def find_next_lunar_solar_birthday_match():
    """
    Finds the next year when the Lunar birthday matches the Solar birthday month and day.
    """
    # Define the person's birth date in the Solar calendar
    birth_year = 1980
    birth_month = 10
    birth_day = 1

    # Convert the Solar birth date to its Lunar equivalent to find the target Lunar birthday
    solar_birth_date = SolarDate(birth_year, birth_month, birth_day)
    lunar_birth_date = solar_birth_date.to_lunar()
    target_lunar_month = lunar_birth_date.month
    target_lunar_day = lunar_birth_date.day

    # Start searching from the year after the birth year
    current_year = birth_year + 1

    while True:
        # Create a SolarDate object for October 1st of the current year
        current_solar_date = SolarDate(current_year, birth_month, birth_day)
        
        # Convert it to the Lunar date for that year
        current_lunar_date = current_solar_date.to_lunar()
        
        # Check if the lunar month and day match the person's lunar birthday
        if (current_lunar_date.month == target_lunar_month and
            current_lunar_date.day == target_lunar_day):
            
            # If they match, we've found the year. Print it and exit the loop.
            print(current_year)
            break
        
        # If not, move to the next year
        current_year += 1

if __name__ == "__main__":
    find_next_lunar_solar_birthday_match()