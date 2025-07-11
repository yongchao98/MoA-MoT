# You may need to install the lunardate library first:
# pip install lunardate

from lunardate import SolarDate

def find_next_matching_birthday():
    """
    Finds the next year after a given birth year where the solar MM-DD
    matches the lunar MM-DD.
    """
    # The person's Solar birthday details
    birth_year = 1980
    target_month = 10
    target_day = 1

    # We start searching from the year after the birth year
    start_year = birth_year + 1

    print(f"Person's Solar Birthday: {birth_year}-{target_month}-{target_day}")
    print(f"Searching for the next year 'Y' where Solar date 'Y-{target_month}-{target_day}' is the same as the Lunar date.\n")

    # Loop through years to find a match. A 100-year range is sufficient.
    for year in range(start_year, start_year + 100):
        # Create a SolarDate object for October 1st of the current year
        solar_date_to_check = SolarDate(year, target_month, target_day)

        # Convert it to its Lunar date equivalent
        lunar_date_equivalent = solar_date_to_check.to_lunar()

        # Check if the lunar month and day match the solar month and day
        if lunar_date_equivalent.month == target_month and lunar_date_equivalent.day == target_day:
            
            # Found the matching year, print the final result
            print(f"Match found in the year: {year}")
            
            print("\nFinal Equation:")
            print(f"Solar Date (Year: {year}, Month: {target_month}, Day: {target_day})")
            print("is equivalent to")
            print(f"Lunar Date (Year: {lunar_date_equivalent.year}, Month: {lunar_date_equivalent.month}, Day: {lunar_date_equivalent.day})")
            
            print(f"\nThe next year is {year}.")
            
            # Return the year and exit the function
            return year

    # In case no match is found in the given range
    print("No matching year found in the next 100 years.")
    return None

# Run the function to get the answer
found_year = find_next_matching_birthday()

# This is a placeholder for the final answer format required by the prompt
# The actual value is printed by the function above.
if found_year:
    # The submission system will extract the answer from this special formatted string.
    # print(f'<<<{found_year}>>>') # This line is for internal processing.
    pass