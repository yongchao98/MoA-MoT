# This script requires the 'lunardate' library.
# You can install it by running: pip install lunardate

from lunardate import LunarDate

def find_next_birthday_coincidence():
    """
    Finds the next year after 1980 when the solar birthday (Oct 1st)
    coincides with the original lunar birthday.
    """
    # The person's solar birth date information.
    birth_solar_year = 1980
    birth_solar_month = 10
    birth_solar_day = 1

    # Step 1: Find the corresponding Lunar birthday for the given solar birth date.
    # This determines the person's fixed lunar birthday (month and day).
    original_lunar_birthday = LunarDate.fromSolarDate(birth_solar_year, birth_solar_month, birth_solar_day)
    target_lunar_month = original_lunar_birthday.month
    target_lunar_day = original_lunar_birthday.day

    # Step 2: Start searching from the year after the birth year.
    year_to_check = birth_solar_year + 1

    while True:
        # Step 3: On each subsequent year's solar birthday (October 1st),
        # convert it to the lunar date for that year.
        current_lunar_date = LunarDate.fromSolarDate(year_to_check, birth_solar_month, birth_solar_day)

        # Step 4: Compare if the calculated lunar date matches the original lunar birthday.
        if current_lunar_date.month == target_lunar_month and current_lunar_date.day == target_lunar_day:
            # If they match, we've found the year.
            print(year_to_check)
            break

        # Move to the next year if no match is found.
        year_to_check += 1

# Execute the function to find and print the answer.
find_next_birthday_coincidence()