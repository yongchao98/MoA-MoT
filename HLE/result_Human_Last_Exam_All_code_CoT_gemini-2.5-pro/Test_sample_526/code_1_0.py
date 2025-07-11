# First, ensure you have the required library installed:
# pip install lunardate

from lunardate import LunarDate, SolarDate

# The person's solar birth date
solar_birth_year = 1980
solar_birth_month = 10
solar_birth_day = 1

# 1. Find the corresponding Lunar birthday for 1980-10-01
try:
    lunar_birth_date = LunarDate.fromSolarDate(solar_birth_year, solar_birth_month, solar_birth_day)
    lunar_birth_month = lunar_birth_date.month
    lunar_birth_day = lunar_birth_date.day

    # 2. Start searching from the year after birth
    year_to_check = solar_birth_year + 1

    while True:
        try:
            # 3. For the current year, find the solar date of the person's lunar birthday
            # This creates a LunarDate object for the lunar birthday in the year we are checking
            lunar_date_in_current_year = LunarDate(year_to_check, lunar_birth_month, lunar_birth_day)

            # Convert it back to its solar equivalent
            solar_equivalent_date = lunar_date_in_current_year.toSolarDate()

            # 4. Check if the solar equivalent matches the original solar birthday (Month and Day)
            if solar_equivalent_date.month == solar_birth_month and solar_equivalent_date.day == solar_birth_day:
                # 5. If it matches, we found the year. Print it and stop.
                print(year_to_check)
                break

        except ValueError:
            # This handles cases where a lunar date might not exist in a specific year
            # (e.g., day 30 in a 29-day lunar month, or issues with leap months).
            # We simply skip to the next year.
            pass

        year_to_check += 1

except ImportError:
    print("The 'lunardate' library is not installed.")
    print("Please install it by running: pip install lunardate")
