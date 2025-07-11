# To run this code, you first need to install the 'lunardate' library.
# You can do this by running the following command in your terminal or command prompt:
# pip install lunardate

import lunardate

def find_next_matching_birthday():
    """
    Finds the next year after 1980 where the Solar date October 1st
    is the same as the Lunar date October 1st.
    """
    # The person was born in 1980.
    birth_year = 1980
    
    # The target birthday month and day are October 1st.
    target_month = 10
    target_day = 1
    
    # We start our search from the year after the birth year.
    start_year = birth_year + 1
    
    # Search for a match in the next 100 years.
    for year in range(start_year, start_year + 100):
        # Create a SolarDate object for October 1st of the current year.
        solar_date = lunardate.SolarDate(year, target_month, target_day)
        
        # Convert it to the corresponding LunarDate.
        lunar_date = solar_date.to_lunar()
        
        # Check if the lunar month and day match the target.
        # We also ensure it's not a leap month by checking 'not lunar_date.isleap'.
        if lunar_date.month == target_month and lunar_date.day == target_day and not lunar_date.isleap:
            # If a match is found, print the year and exit the function.
            print(year)
            return

# Run the function to find and print the result.
find_next_matching_birthday()
