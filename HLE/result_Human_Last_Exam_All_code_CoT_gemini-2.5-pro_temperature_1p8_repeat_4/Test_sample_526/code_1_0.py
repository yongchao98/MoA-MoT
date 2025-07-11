# The user needs to install the 'lunardate' library first.
# You can install it by running: pip install lunardate
from lunardate import LunarDate, SolarDate

def find_next_matching_birthday_year():
    """
    Finds the next year after 1980 where the Lunar birthday for a person
    born on 1980-10-01 falls on the same Solar date (10-01).
    """
    # 1. Define the birth date in the Solar calendar.
    birth_solar_year = 1980
    birth_solar_month = 10
    birth_solar_day = 1

    # 2. Convert the Solar birth date to its Lunar equivalent to find the Lunar month and day.
    try:
        birth_lunar_date = LunarDate.fromSolarDate(birth_solar_year, birth_solar_month, birth_solar_day)
        birth_lunar_month = birth_lunar_date.month
        birth_lunar_day = birth_lunar_date.day
    except ImportError:
        print("The 'lunardate' library is not installed.")
        print("Please install it using: pip install lunardate")
        return
    except Exception as e:
        print(f"An error occurred: {e}")
        return

    # 3. Loop through the years starting from the year after birth.
    # We set a reasonable upper limit to prevent an infinite loop.
    for year in range(birth_solar_year + 1, birth_solar_year + 100):
        try:
            # 4. For the current year, create a LunarDate object with the person's lunar birth month and day.
            # This represents their Lunar birthday in the current year.
            # The lunardate library handles leap months automatically.
            current_lunar_birthday = LunarDate(year, birth_lunar_month, birth_lunar_day)

            # Convert this Lunar birthday to its Solar date equivalent.
            corresponding_solar_date = current_lunar_birthday.toSolarDate()

            # 5. Check if the resulting Solar date matches the original birth month and day.
            if (corresponding_solar_date.month == birth_solar_month and
                corresponding_solar_date.day == birth_solar_day):
                # 6. If it matches, we have found the answer.
                print(corresponding_solar_date.year)
                return

        except ValueError:
            # This handles cases where a lunar month/day doesn't exist in a particular year (e.g., leap months).
            continue

if __name__ == '__main__':
    find_next_matching_birthday_year()
<<<1999>>>