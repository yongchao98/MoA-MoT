import sys

# Attempt to import the required library for Lunar-Solar calendar conversions.
try:
    from lunardate import Solar, Lunar
except ImportError:
    # If the library is not found, print an error message and exit.
    print("Error: The 'lunardate' library is required to run this script.", file=sys.stderr)
    print("Please install it using pip: pip install lunardate", file=sys.stderr)
    sys.exit(1)

def find_next_birthday_match():
    """
    Finds the next year after 1980 where the Lunar birthday (for a person
    born on 1980-10-01) falls on the same Solar date (10-01).
    """
    # The person's birth date in the Solar calendar.
    birth_year_solar = 1980
    birth_month_solar = 10
    birth_day_solar = 1

    # Convert the solar birth date to the lunar calendar to find the lunar birthday.
    # Solar(year, month, day)
    solar_birth_date = Solar(birth_year_solar, birth_month_solar, birth_day_solar)
    lunar_birth_date = solar_birth_date.to_lunar()
    
    # These are the lunar month and day we are looking for in subsequent years.
    # For 1980-10-01, this corresponds to Lunar Month 8, Day 23.
    target_lunar_month = lunar_birth_date.month
    target_lunar_day = lunar_birth_date.day

    # Start searching from the year after birth.
    current_year = birth_year_solar + 1

    # Loop indefinitely until we find the matching year.
    while True:
        try:
            # Create a lunar date for the current year using the target lunar month and day.
            # Lunar(year, month, day)
            lunar_date_in_current_year = Lunar(current_year, target_lunar_month, target_lunar_day)
            
            # Convert this lunar date back to its solar equivalent.
            solar_equivalent_date = lunar_date_in_current_year.to_solar()

            # Check if the resulting solar date's month and day match the original solar birthday.
            if solar_equivalent_date.month == birth_month_solar and solar_equivalent_date.day == birth_day_solar:
                # We have found the year. Print the final answer.
                # The final equation is essentially finding Y such that:
                # Solar(Lunar(Y, 8, 23)) == YYYY-10-01
                # Where Lunar(1980, 8, 23) was the original birthday.
                print(f"Birth Date (Solar): {birth_year_solar}-{birth_month_solar}-{birth_day_solar}")
                print(f"Lunar Birthday: Month {target_lunar_month}, Day {target_lunar_day}")
                print(f"Found a match in year: {current_year}")
                print(f"In {current_year}, Lunar {target_lunar_month}-{target_lunar_day} is Solar {solar_equivalent_date.year}-{solar_equivalent_date.month}-{solar_equivalent_date.day}")
                print("\nFinal Answer:")
                print(current_year)
                return current_year

        except ValueError:
            # This exception occurs if the lunar date does not exist in the current year
            # (e.g., a leap month date in a non-leap year). We simply continue to the next year.
            pass
        
        # If no match was found, move to the next year.
        current_year += 1

        # A safeguard to prevent an accidental infinite loop.
        if current_year > birth_year_solar + 150:
            print("Failed to find a match within 150 years.", file=sys.stderr)
            break

# Execute the function.
if __name__ == "__main__":
    find_next_birthday_match()