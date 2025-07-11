# The user needs to install the 'lunardate' package first.
# They can do this by running: pip install lunardate

from lunardate import SolarDate

def find_next_matching_birthday():
    """
    Finds the next year after 1980 where the Solar date 10-01
    is also the Lunar date 10-01.
    """
    birth_year = 1980
    target_month = 10
    target_day = 1

    # Start searching from the year after birth
    for year in range(birth_year + 1, birth_year + 150):
        # Create a solar date object for October 1st of the current year
        solar_date = SolarDate(year, target_month, target_day)
        
        # Convert the solar date to its lunar equivalent
        lunar_date = solar_date.to_lunar()
        
        # Check if the lunar month and day match the target month and day
        if lunar_date.month == target_month and lunar_date.day == target_day:
            # If they match, we have found the year.
            # Print the final equation/result details
            print(f"The person was born on Solar date 1980-{target_month:02d}-{target_day:02d}.")
            print(f"The goal is to find the next year where Solar 'MM-DD' is the same as Lunar 'MM-DD' (i.e., {target_month:02d}-{target_day:02d}).")
            print(f"\nFound a match in the year: {year}")
            print(f"In this year, the Solar date {year}-{solar_date.month:02d}-{solar_date.day:02d} corresponds to the Lunar date {lunar_date.year}-{lunar_date.month:02d}-{lunar_date.day:02d}.")
            print("\nThe final answer is the year when this occurs.")
            print(f"{year}")
            return year

    # If no year is found in the range
    print(f"No year found within the search range where Solar {target_month}-{target_day} matches Lunar {target_month}-{target_day}.")
    return None

# Run the function to find and print the answer.
found_year = find_next_matching_birthday()

# The final answer in the required format will be derived from the printed output.
# The code above prints the answer. For the final submission format:
# print(f"<<<{found_year}>>>")