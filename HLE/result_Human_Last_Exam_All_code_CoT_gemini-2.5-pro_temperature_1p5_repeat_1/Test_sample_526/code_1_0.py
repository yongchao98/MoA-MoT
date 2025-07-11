# You may need to install the required library first:
# pip install pysxtwl

import sxtwl

def find_next_matching_birthday():
    """
    Finds the next year after a given birth year where the Solar and Lunar birthdays match.
    """
    solar_birth_year = 1980
    solar_birth_month = 10
    solar_birth_day = 1

    # Get the time zone utility for Chinese calendar calculations
    lunar = sxtwl.Lunar()

    # Convert the solar birth date to its lunar equivalent
    try:
        birth_day_info = lunar.getDayBySolar(solar_birth_year, solar_birth_month, solar_birth_day)
    except Exception as e:
        print(f"Error converting the initial birth date: {e}")
        return

    # Store the specific lunar month and day from the birth date
    lunar_birth_month = birth_day_info.Lmc
    lunar_birth_day = birth_day_info.Ldc
    # Check if the birth month was a leap month
    is_birth_month_leap = birth_day_info.Lleap

    # Start searching from the year after birth
    current_year = solar_birth_year + 1
    # We set a reasonable search limit, e.g., 100 years
    search_limit = solar_birth_year + 100

    while current_year < search_limit:
        try:
            # For the current year, convert the lunar birthday back to a solar date
            solar_date_candidate = lunar.getDayByLunar(current_year, lunar_birth_month, lunar_birth_day, is_birth_month_leap)

            # Check if the resulting solar date is October 1st
            if solar_date_candidate.Smc == solar_birth_month and solar_date_candidate.Sdc == solar_birth_day:
                # If it matches, we have found the year
                print(current_year)
                return
        except Exception:
            # This exception can occur if the lunar date is invalid for the current year
            # (e.g., a leap month that doesn't exist, or day 30 in a 29-day month).
            # We simply ignore this year and continue to the next one.
            pass

        current_year += 1

    print("Could not find a matching year within the search limit.")

if __name__ == '__main__':
    find_next_matching_birthday()