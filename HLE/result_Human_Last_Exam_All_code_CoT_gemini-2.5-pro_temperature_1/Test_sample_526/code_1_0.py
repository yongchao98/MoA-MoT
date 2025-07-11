# The user may need to install the required library first.
# You can do this by running: pip install lunardate

from lunardate import SolarDate, LunarDate

def find_next_matching_birthday():
    """
    Finds the next year after 1980 where the Solar date 10-01
    matches the Lunar date 10-01.
    """
    birth_year = 1980
    target_month = 10
    target_day = 1

    # Start checking from the year after birth
    for year in range(birth_year + 1, birth_year + 100):
        # Create a solar date object for Oct 1st of the current year
        solar_date = SolarDate(year, target_month, target_day)
        
        # Convert the solar date to its lunar equivalent
        lunar_date = solar_date.to_lunar()
        
        # Check if the lunar month and day match the target
        if lunar_date.month == target_month and lunar_date.day == target_day:
            print(f"The person was born on 1980-10-01.")
            print(f"The target is to find a year > 1980 where the Solar date {target_month:02d}-{target_day:02d} is also the Lunar date {target_month:02d}-{target_day:02d}.")
            print(f"In the year {year}, the Solar date {year}-{target_month:02d}-{target_day:02d} corresponds to the Lunar date {lunar_date.year}-{lunar_date.month:02d}-{lunar_date.day:02d}.")
            print(f"Therefore, the next year is: {year}")
            return year

    print("No matching year found in the next 100 years.")
    return None

if __name__ == '__main__':
    result_year = find_next_matching_birthday()
    if result_year:
        print(f"<<<{result_year}>>>")
