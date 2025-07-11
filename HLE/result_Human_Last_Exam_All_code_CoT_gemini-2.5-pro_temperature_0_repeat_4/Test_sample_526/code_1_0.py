# The 'lunardate' library is required for this script.
# If you don't have it, you can install it by running: pip install lunardate
try:
    from lunardate import LunarDate
except ImportError:
    print("Error: The 'lunardate' library is not installed.")
    print("Please install it by running the command: pip install lunardate")
    exit()

# The person's solar birth date details
birth_year = 1980
birth_month = 10
birth_day = 1

# 1. Find the corresponding lunar date for the birth date
birth_lunar_date = LunarDate.fromSolarDate(birth_year, birth_month, birth_day)
birth_lunar_month = birth_lunar_date.month
birth_lunar_day = birth_lunar_date.day

# 2. Start checking from the year after the birth year
year_to_check = birth_year + 1

# Loop until we find the matching year
while True:
    # 3. Get the lunar date for October 1st of the current year being checked
    current_lunar_date = LunarDate.fromSolarDate(year_to_check, birth_month, birth_day)

    # 4. Check if the lunar month and day match the person's lunar birthday
    if current_lunar_date.month == birth_lunar_month and current_lunar_date.day == birth_lunar_day:
        # 5. We found the year. Print the details of the match.
        print(f"The person's solar birthday is {birth_year}-{birth_month:02d}-{birth_day:02d}.")
        print(f"The corresponding lunar birthday is Month {birth_lunar_month}, Day {birth_lunar_day}.")
        print(f"The next year when the solar date {year_to_check}-{birth_month:02d}-{birth_day:02d} matches this lunar birthday is:")
        print(year_to_check)
        break

    # Move to the next year
    year_to_check += 1

    # Add a safeguard to prevent an infinite loop, though a match is expected.
    if year_to_check > birth_year + 200:
        print("Could not find a matching year within 200 years.")
        break