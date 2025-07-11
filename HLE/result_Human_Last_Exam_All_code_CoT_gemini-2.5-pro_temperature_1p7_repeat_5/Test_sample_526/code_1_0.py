# First, you may need to install the lunardate library:
# pip install lunardate

import lunardate
from datetime import date

def find_next_birthday_match():
    """
    Finds the next year when a person's Lunar and Solar birthdays (MM-DD) align.
    """
    # The person's Solar birth date
    birth_year = 1980
    birth_month = 10
    birth_day = 1

    # Step 1: Find the Lunar birthday (month and day) for the given Solar birth date.
    try:
        solar_birth_date = date(birth_year, birth_month, birth_day)
        lunar_birth_date = lunardate.LunarDate.fromSolarDate(
            solar_birth_date.year,
            solar_birth_date.month,
            solar_birth_date.day
        )
        target_lunar_month = lunar_birth_date.month
        target_lunar_day = lunar_birth_date.day
    except Exception as e:
        print(f"Could not determine the Lunar birth date from the provided Solar date. Error: {e}")
        return

    # Step 2: Loop through subsequent years to find the next match.
    # We start searching from the year after birth.
    for year_to_check in range(birth_year + 1, birth_year + 100):
        try:
            # Step 3: For the current year, determine the Solar date of the Lunar birthday.
            lunar_date_in_new_year = lunardate.LunarDate(year_to_check, target_lunar_month, target_lunar_day)
            corresponding_solar_date = lunar_date_in_new_year.toSolarDate()

            # Step 4: Check if the corresponding Solar date's month and day match the original.
            if corresponding_solar_date.month == birth_month and corresponding_solar_date.day == birth_day:
                # A match is found. This is the year we are looking for.
                found_year = corresponding_solar_date.year
                print(found_year)
                return

        except (ValueError, KeyError):
            # This handles cases like leap lunar months or invalid dates for a specific year.
            # We can safely skip these years and continue the search.
            continue
    
    print("A matching year could not be found within the next 100 years.")

if __name__ == '__main__':
    find_next_birthday_match()