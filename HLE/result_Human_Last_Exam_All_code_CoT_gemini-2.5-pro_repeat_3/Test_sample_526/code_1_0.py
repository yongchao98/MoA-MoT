import datetime
from lunardate import LunarDate

def find_next_matching_birthday():
    """
    Finds the next year after a given birth year where the Solar birthday
    and the Lunar birthday fall on the same month and day again.
    """
    # Define the person's solar birth date
    birth_year = 1980
    birth_month = 10
    birth_day = 1

    # 1. Find the original Lunar birthday (MM-DD)
    try:
        original_lunar_date = LunarDate.fromSolarDate(birth_year, birth_month, birth_day)
        target_lunar_month = original_lunar_date.month
        target_lunar_day = original_lunar_date.day
    except ImportError:
        print("Please install the 'lunardate' library by running: pip install lunardate")
        return

    # 2. Iterate through years starting from the year after birth
    # A reasonable upper limit is 100 years.
    for year in range(birth_year + 1, birth_year + 101):
        # Get the lunar date for the same solar day (Oct 1st) in the current year
        current_solar_date = datetime.date(year, birth_month, birth_day)
        current_lunar_date = LunarDate.fromSolarDate(current_solar_date.year, current_solar_date.month, current_solar_date.day)

        # 3. Check if the lunar month and day match the original lunar birthday
        if current_lunar_date.month == target_lunar_month and current_lunar_date.day == target_lunar_day:
            # Found the matching year
            print(f"The person was born on {birth_year}-{birth_month:02d}-{birth_day:02d} (Solar).")
            print(f"Their Lunar birthday is Month {target_lunar_month}, Day {target_lunar_day}.")
            print(f"The next year their Solar birthday ({year}-{birth_month:02d}-{birth_day:02d}) falls on the same Lunar date is:")
            print(year)
            return

    print("No matching year found within 100 years.")

if __name__ == '__main__':
    find_next_matching_birthday()
