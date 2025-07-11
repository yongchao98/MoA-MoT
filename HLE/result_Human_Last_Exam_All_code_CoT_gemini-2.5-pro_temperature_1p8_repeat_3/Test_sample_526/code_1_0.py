import lunardate
from lunardate import LunarDate

def find_next_birthday_match():
    """
    Finds the next year after 1980 where the Lunar and Solar birthdays
    (October 1st) fall on the same day.
    """
    # Step 1: Define the solar birth date
    birth_year = 1980
    birth_month = 10
    birth_day = 1

    # Step 2: Convert the solar birth date to find the lunar birthday (month and day)
    lunar_birth_date = LunarDate.fromSolarDate(birth_year, birth_month, birth_day)
    target_lunar_month = lunar_birth_date.month
    target_lunar_day = lunar_birth_date.day

    print(f"The solar birth date is: {birth_year}-{birth_month}-{birth_day}")
    print(f"The corresponding lunar birthday is: Month {target_lunar_month}, Day {target_lunar_day}\n")

    # Step 3: Iterate through years starting from the year after birth
    start_year = birth_year + 1
    # Set a reasonable upper limit for the search
    end_year = start_year + 100

    print("Searching for the next year where the two birthdays align...")

    for year in range(start_year, end_year):
        try:
            # Construct the lunar date for the current year
            current_lunar_date = LunarDate(year, target_lunar_month, target_lunar_day)
            
            # Convert it back to a solar date
            corresponding_solar_date = current_lunar_date.toSolarDate()

            # Step 4: Check if the solar month and day match the original birth date
            if corresponding_solar_date.month == birth_month and corresponding_solar_date.day == birth_day:
                print(f"Found a match in the year {year}.")
                print(f"In {year}, Lunar date {target_lunar_month}-{target_lunar_day} corresponds to Solar date {corresponding_solar_date.month}-{corresponding_solar_date.day}.")
                print("\nThe next year is:")
                print(year)
                return

        except ValueError:
            # This handles cases where a lunar date doesn't exist in a particular year,
            # e.g., the 30th day of a lunar month with only 29 days.
            continue

if __name__ == '__main__':
    find_next_birthday_match()