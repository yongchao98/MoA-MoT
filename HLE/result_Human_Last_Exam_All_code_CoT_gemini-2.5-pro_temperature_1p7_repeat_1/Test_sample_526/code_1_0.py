# To run this code, you first need to install the `lunardate` library.
# You can do this by running the following command in your terminal or command prompt:
# pip install lunardate

import lunardate

def find_next_matching_lunar_birthday():
    """
    A person is born on 1980-10-01. This function finds the next year
    where the Solar calendar date of October 1st is the same as the
    Lunar calendar date (i.e., the 1st day of the 10th lunar month).
    """
    # The person's birth year
    birth_year = 1980
    
    # The target Solar month and day
    target_month = 10
    target_day = 1

    # We start searching from the year after the birth year
    year_to_check = birth_year + 1

    # Loop indefinitely until we find the matching year
    while True:
        # Convert the Solar date (YYYY-10-01) to the Lunar date for the current year
        lunar_date = lunardate.LunarDate.fromSolarDate(year_to_check, target_month, target_day)

        # Check if the converted Lunar date is also October 1st
        if lunar_date.month == target_month and lunar_date.day == target_day:
            # We found the matching year. Print the result.
            print(f"The person was born in the year {birth_year}.")
            print(f"The target solar and lunar date is month: {target_month}, day: {target_day}.")
            print(f"The next year where this match occurs is: {year_to_check}")
            return

        # If it's not a match, move on to the next year
        year_to_check += 1

# --- Main execution ---
# Run the function to find and print the answer
find_next_matching_lunar_birthday()