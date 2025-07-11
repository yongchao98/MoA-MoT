import sys
try:
    from lunardate import LunarDate, SolarDate
except ImportError:
    print("The 'lunardate' library is required. Please install it by running: pip install lunardate", file=sys.stderr)
    sys.exit(1)

# The person's solar birthday
birth_year = 1980
birth_month = 10
birth_day = 1

# 1. Convert the solar birth date to a lunar date to find the lunar birthday
solar_birth_date = SolarDate(birth_year, birth_month, birth_day)
lunar_birth_date = solar_birth_date.to_lunar()
target_lunar_month = lunar_birth_date.month
target_lunar_day = lunar_birth_date.day

# 2. Iterate through years starting from the year after birth
start_year = birth_year + 1
# Set a reasonable upper limit for the search
end_year = start_year + 100 

for year_to_check in range(start_year, end_year):
    try:
        # Create a LunarDate object for the target lunar birthday in the current year
        lunar_date_in_year = LunarDate(year_to_check, target_lunar_month, target_lunar_day)

        # Convert this lunar date back to a solar date
        corresponding_solar_date = lunar_date_in_year.to_solar()

        # 3. Check if the resulting solar date matches the original solar birthday's month and day
        if corresponding_solar_date.month == birth_month and corresponding_solar_date.day == birth_day:
            # 4. If it matches, print the year and exit the loop
            print(year_to_check)
            break
    except ValueError:
        # Some years might not have the target lunar date (e.g., leap months can cause this).
        # We can safely ignore these years and continue the search.
        continue
<<<1999>>>