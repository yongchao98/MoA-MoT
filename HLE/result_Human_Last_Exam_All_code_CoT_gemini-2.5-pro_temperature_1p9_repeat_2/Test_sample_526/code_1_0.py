import lunardate
import datetime

def find_next_matching_birthday_year():
    """
    Finds the next year when the Lunar birthday falls on the same Solar month and day
    as the original Solar birthday.
    """
    # Define the person's Solar birth date information
    birth_year = 1980
    birth_month = 10
    birth_day = 1

    # 1. Convert the Solar birth date to its Lunar equivalent to find the fixed Lunar birthday.
    try:
        solar_birth_date = lunardate.Solar(birth_year, birth_month, birth_day)
        lunar_birth_date = solar_birth_date.toLunar()
    except (ValueError, AttributeError) as e:
        print(f"Error converting the birth date: {e}")
        return

    lunar_birth_month = lunar_birth_date.month
    lunar_birth_day = lunar_birth_date.day
    is_leap_month_birth = lunar_birth_date.isleap

    # 2. Start searching from the year after birth.
    year_to_check = birth_year + 1

    while True:
        try:
            # 3. Create the Lunar date for the person's birthday in the current year of the loop.
            # We must use the 'isleap' flag from their birth year to handle birthdays in leap months correctly.
            current_lunar_birthday = lunardate.Lunar(year_to_check, lunar_birth_month, lunar_birth_day, is_leap_month_birth)
            
            # Convert this Lunar date back to its Solar equivalent for that year.
            corresponding_solar_date = current_lunar_birthday.toSolar()

            # 4. Check if the month and day of the new Solar date match the original birth month and day.
            if corresponding_solar_date.month == birth_month and corresponding_solar_date.day == birth_day:
                # We found the matching year.
                print(f"The next year after {birth_year} where the Lunar birthday and Solar birthday match is: {year_to_check}")
                break # Exit the loop once the answer is found.

        except ValueError:
            # This handles years where the specific Lunar date doesn't exist.
            # For example, a leap month birthday in a year with no such leap month,
            # or day 30 in a lunar month that only has 29 days.
            pass
            
        year_to_check += 1
        
        # Add a safeguard to prevent an infinite loop in case no match is found within a reasonable time.
        if year_to_check > birth_year + 100:
            print("Could not find a matching year within 100 years.")
            break

find_next_matching_birthday_year()