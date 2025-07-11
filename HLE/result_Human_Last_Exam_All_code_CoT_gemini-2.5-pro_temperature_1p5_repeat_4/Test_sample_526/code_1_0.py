import lunardate
from lunardate import Lunar, Solar

def find_next_matching_birthday():
    """
    Finds the next year after 1980 where the Lunar birthday corresponding to
    1980-10-01 falls on the same Solar date (10-01).
    """
    # Define the Solar birth date
    birth_year = 1980
    birth_month = 10
    birth_day = 1

    # 1. Convert the Solar birth date to its Lunar equivalent
    try:
        solar_birth_date = Solar(birth_year, birth_month, birth_day)
        lunar_birth_date = solar_birth_date.to_lunar()
        target_lunar_month = lunar_birth_date.month
        target_lunar_day = lunar_birth_date.day
    except ImportError:
        print("The 'lunardate' library is required. Please install it using 'pip install lunardate'")
        return

    # 2. Iterate through subsequent years to find a match
    # A search window of 100 years is sufficient.
    for year_to_check in range(birth_year + 1, birth_year + 101):
        try:
            # 3. Create a Lunar date object for the current year being checked
            lunar_date_in_year = Lunar(year_to_check, target_lunar_month, target_lunar_day)

            # Convert this Lunar date back to its Solar equivalent
            solar_date_in_year = lunar_date_in_year.to_solar()

            # 4. Check if the resulting Solar date's month and day match the original birth date
            if solar_date_in_year.month == birth_month and solar_date_in_year.day == birth_day:
                
                # 5. Found the answer. Print the final equation-like result.
                print(f"Solar Birthday: {birth_year}-{birth_month}-{birth_day}")
                print(f"Corresponding Lunar Birthday: Month {target_lunar_month}, Day {target_lunar_day}")
                print("-" * 30)
                print(f"Found match in year: {year_to_check}")
                print(f"Lunar Date ({year_to_check}, {target_lunar_month}, {target_lunar_day}) -> Solar Date ({solar_date_in_year.year}, {solar_date_in_year.month}, {solar_date_in_year.day})")
                print(f"\nThe next year is: {year_to_check}")
                return year_to_check
                
        except ValueError:
            # This handles cases where a lunar date is invalid for a given year (e.g., leap months)
            continue
    
    print("No matching year was found in the next 100 years.")
    return None

# Execute the function to find and print the answer
found_year = find_next_matching_birthday()

# Final answer block as requested by the format
if found_year:
    print(f"<<<{found_year}>>>")
