# First, you might need to install the required library:
# pip install lunardate

from lunardate import SolarDate

# The person's Solar birthday is October 1st.
target_month = 10
target_day = 1

# The person was born in 1980. We need to find the *next* year.
start_year = 1981

# We will loop through the years to find a match.
# A reasonable upper limit is set to prevent an infinite loop.
for year in range(start_year, start_year + 100):
    # Create a SolarDate object for October 1st of the current year.
    solar_date = SolarDate(year, target_month, target_day)
    
    # Convert the Solar date to its Lunar equivalent.
    lunar_date = solar_date.to_lunar()
    
    # Check if the Lunar month and day match the target birthday.
    # We also need to ensure it's not a leap month (is_leap_month is False).
    if lunar_date.month == target_month and lunar_date.day == target_day and not lunar_date.is_leap_month:
        # If they match, we have found the year.
        print(f"The person was born on 1980-10-01.")
        print(f"The next year when their Solar and Lunar birthdays (MM-DD) align is:")
        print(year)
        break