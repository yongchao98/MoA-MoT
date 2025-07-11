# You may need to install the 'lunardate' library first by running:
# pip install lunardate

from lunardate import LunarDate

# 1. Define the solar birth date details
birth_year = 1980
birth_month = 10
birth_day = 1

# 2. Find the person's Lunar birthday (Month and Day)
birth_lunar_date = LunarDate.fromSolarDate(birth_year, birth_month, birth_day)
target_lunar_month = birth_lunar_date.month
target_lunar_day = birth_lunar_date.day

print(f"The person's solar birthday is {birth_year}-{birth_month:02d}-{birth_day:02d}.")
print(f"This corresponds to the lunar birthday: Month = {target_lunar_month}, Day = {target_lunar_day}.")
print("-" * 40)
print("Searching for the next year with the same solar and lunar birthday (MM-DD)...")
print("-" * 40)


# 3. Iterate through subsequent years to find the next match
year_to_check = birth_year + 1
while True:
    # Get the lunar date for October 1st of the year we are checking
    current_lunar_date = LunarDate.fromSolarDate(year_to_check, birth_month, birth_day)

    # Check if the lunar month and day match the target birthday
    if current_lunar_date.month == target_lunar_month and current_lunar_date.day == target_lunar_day:
        found_year = year_to_check
        
        print(f"MATCH FOUND!")
        # The prompt requires outputting the numbers in the final equation.
        # Here we show the final matching equation/check.
        print(f"Year: {found_year}")
        print(f"Solar Date: {found_year}-{birth_month}-{birth_day}")
        print(f"Lunar Date: Month={current_lunar_date.month}, Day={current_lunar_date.day}")
        print(f"This matches the birth lunar date: Month={target_lunar_month}, Day={target_lunar_day}")
        
        print("\nFinal Answer (Year):")
        print(found_year)
        break

    # Add a safety break to prevent an infinite loop, though a match is expected.
    if year_to_check > birth_year + 100:
        print("No match found within 100 years.")
        break
        
    year_to_check += 1