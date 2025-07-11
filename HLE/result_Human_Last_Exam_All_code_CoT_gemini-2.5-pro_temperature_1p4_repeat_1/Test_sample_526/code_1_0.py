import sys

# The 'lunardate' library is required for this script.
# If you don't have it, please install it first by running:
# pip install lunardate

try:
    from lunardate import SolarDate, LunarDate
except ImportError:
    print("Error: The 'lunardate' library is required but not found.")
    print("Please install it using the command: pip install lunardate")
    sys.exit(1)

# The person's solar birth date
birth_year_solar = 1980
birth_month_solar = 10
birth_day_solar = 1

# Step 1: Convert the solar birth date to its lunar equivalent
# to find the target lunar month and day.
solar_birth_date = SolarDate(birth_year_solar, birth_month_solar, birth_day_solar)
lunar_birth_date = solar_birth_date.toLunarDate()
target_lunar_month = lunar_birth_date.month
target_lunar_day = lunar_birth_date.day

# Step 2: Start searching from the year after birth
year_to_check = birth_year_solar + 1

while True:
    try:
        # Step 3: For the current year, create a lunar date using the target month and day
        lunar_date_in_year = LunarDate(year_to_check, target_lunar_month, target_lunar_day)
        
        # Convert this lunar date back to its solar equivalent
        solar_equivalent_date = lunar_date_in_year.toSolarDate()
        
        # Step 4: Check if the resulting solar date matches the original birth month and day
        if solar_equivalent_date.month == birth_month_solar and solar_equivalent_date.day == birth_day_solar:
            # If a match is found, print the year and stop searching
            print(year_to_check)
            break
            
    except (ValueError, IndexError):
        # This error can occur if the lunar date is invalid for a specific year
        # (e.g., trying to find day 30 in a 29-day lunar month).
        # We can safely ignore it and proceed to the next year.
        pass

    # Move to the next year
    year_to_check += 1