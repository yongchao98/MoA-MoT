import sys
# To run this code, you may need to install the lunardate library first:
# pip install lunardate
try:
    from lunardate import LunarDate
except ImportError:
    print("Error: The 'lunardate' library is not installed.", file=sys.stderr)
    print("Please install it by running: pip install lunardate", file=sys.stderr)
    sys.exit(1)

# 1. Define the Solar birth date.
birth_year_solar = 1980
birth_month_solar = 10
birth_day_solar = 1

# 2. Convert the Solar date to its Lunar equivalent to find the Lunar birthday.
lunar_birth_date = LunarDate.fromSolarDate(birth_year_solar, birth_month_solar, birth_day_solar)
target_lunar_month = lunar_birth_date.month
target_lunar_day = lunar_birth_date.day

# 3. Iterate through years starting from the one after birth.
year_to_check = birth_year_solar + 1
found_year = None

# Set a reasonable search limit.
while year_to_check < 2100:
    try:
        # For the current year, find the Solar date of the person's Lunar birthday.
        solar_equivalent_date = LunarDate(year_to_check, target_lunar_month, target_lunar_day).toSolarDate()
        
        # 4. Check if the resulting Solar date matches the original birthday (MM-DD).
        if solar_equivalent_date.month == birth_month_solar and solar_equivalent_date.day == birth_day_solar:
            found_year = year_to_check
            break  # Exit the loop as we've found the first match.
            
    except (ValueError, IndexError):
        # This handles cases where a lunar date is invalid in a given year
        # (e.g., trying to find day 30 in a 29-day lunar month).
        pass
        
    year_to_check += 1

# 5. Print the final result.
if found_year:
    print(found_year)
else:
    print("A matching year could not be found in the specified range.")
