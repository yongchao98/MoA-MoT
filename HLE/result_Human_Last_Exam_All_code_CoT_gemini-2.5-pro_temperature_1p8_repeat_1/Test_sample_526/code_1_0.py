# First, you may need to install the library: pip install lunardate
from lunardate import LunarDate
from datetime import date

def find_next_matching_birthday():
    """
    Finds the next year after 1980 where the Solar date 10-01
    is also the Lunar date 10-01.
    """
    birth_year = 1980
    target_solar_month = 10
    target_solar_day = 1
    
    # We will search from the year after birth up to 100 years later.
    for year in range(birth_year + 1, birth_year + 101):
        # Create a solar date for October 1st of the current year in the loop.
        solar_date = date(year, target_solar_month, target_solar_day)
        
        # Convert the solar date to a lunar date.
        lunar_date = LunarDate.fromSolarDate(year, target_solar_month, target_solar_day)
        
        # Check if the lunar month and day match the target solar month and day.
        if lunar_date.month == target_solar_month and lunar_date.day == target_solar_day:
            print(f"The next year when Solar and Lunar birthdays match is: {year}")
            print(f"In that year, the Solar date {solar_date.month}-{solar_date.day} is the same as the Lunar date {lunar_date.month}-{lunar_date.day}.")
            # The final answer in the required format
            print(f"<<<{year}>>>")
            return

    print("No matching year was found in the next 100 years.")

if __name__ == '__main__':
    find_next_matching_birthday()
