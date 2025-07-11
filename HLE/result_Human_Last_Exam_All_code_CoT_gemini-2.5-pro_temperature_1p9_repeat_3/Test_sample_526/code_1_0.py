# First, you may need to install the lunardate library.
# You can do this by running the following command in your terminal or command prompt:
# pip install lunardate

from lunardate import SolarDate, LunarDate

# The person's birth year
birth_year = 1980

# The target Solar birthday month and day
target_solar_month = 10
target_solar_day = 1

# The target Lunar birthday month and day, which are the same
target_lunar_month = 10
target_lunar_day = 1

# We will start searching from the year after the birth year
start_year = birth_year + 1

# Set a reasonable range for the search (e.g., 100 years)
end_year = start_year + 100

found_year = None

for year in range(start_year, end_year):
    # Create a SolarDate object for October 1st of the current year
    solar_date = SolarDate(year, target_solar_month, target_solar_day)
    
    # Convert the Solar date to its corresponding Lunar date
    lunar_date = solar_date.to_lunar()
    
    # Check if the Lunar month and day match the target
    if lunar_date.month == target_lunar_month and lunar_date.day == target_lunar_day:
        found_year = year
        print(f"Found a match in the year: {found_year}")
        print(f"Solar Date: {solar_date.year}-{solar_date.month}-{solar_date.day}")
        print(f"Lunar Date: {lunar_date.year}-{lunar_date.month}-{lunar_date.day}")
        print(f"\nThe next year since their birth that their Lunar calendar birthday is the same as the Solar Calendar's birthday is: {found_year}")
        break

if found_year is None:
    print(f"No matching year found between {start_year} and {end_year-1}.")

<<<2039>>>