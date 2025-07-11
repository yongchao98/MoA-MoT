# The 'lunardate' library is required to run this script.
# If you don't have it, you can install it by running: pip install lunardate

from lunardate import SolarDate

# The person was born in 1980. We need to find the *next* year.
start_year = 1980
year = start_year + 1

# We are looking for a match where Solar date MM-DD is the same as Lunar date MM-DD.
target_month = 10
target_day = 1

while True:
    # Create a SolarDate object for October 1st of the current year.
    solar_date = SolarDate(year, target_month, target_day)
    
    # Convert the Solar date to its Lunar equivalent.
    lunar_date = solar_date.to_lunar()
    
    # Check if the Lunar date is also October 1st.
    if lunar_date.month == target_month and lunar_date.day == target_day:
        # We found the matching year.
        # Print the numbers that make up the final result year.
        year_str = str(year)
        for digit in year_str:
            print(digit)

        # To be extra clear, we also print the full year on a new line.
        print(year)
        break # Exit the loop once the answer is found.
    
    # If no match, move to the next year.
    year += 1