# First, ensure you have the lunardate library installed:
# pip install lunardate

import datetime
from lunardate import LunarDate

def find_next_matching_lunar_solar_date():
    """
    Finds the next year after 1980 where the Solar date of October 1st
    is also the Lunar date of October 1st.
    """
    # The person was born in 1980, so we start the search from the next year.
    start_year = 1981
    target_month = 10
    target_day = 1

    # We will loop through a reasonable range of years to find the match.
    for year in range(start_year, 2101):
        # Convert the Solar date (YYYY-10-01) to its corresponding Lunar date.
        lunar_date = LunarDate.fromSolarDate(year, target_month, target_day)

        # Check if the Lunar month and day match our target.
        if lunar_date.month == target_month and lunar_date.day == target_day:
            # If they match, we have found the year.
            # We use a variable to make it clear what number we are printing.
            result_year = year
            print(result_year)
            return

    # This part would be executed if no match is found in the given range.
    print("No matching year found before 2101.")

if __name__ == "__main__":
    find_next_matching_lunar_solar_date()
