# First, you may need to install the required library. You can do this by running the following command in your terminal:
# pip install lunardate

from lunardate import SolarDate

def find_next_matching_birthday():
    """
    This function finds the next year after 1980 where the Lunar calendar birthday (MM-DD)
    matches the Solar calendar birthday (10-01).
    """
    # 1. Define the Solar birth date.
    birth_year = 1980
    birth_month = 10
    birth_day = 1
    
    solar_birth_date = SolarDate(birth_year, birth_month, birth_day)

    # 2. Find the corresponding Lunar birthday (month and day) to use as our target.
    lunar_birth_date = solar_birth_date.to_lunar()
    target_lunar_month = lunar_birth_date.month
    target_lunar_day = lunar_birth_date.day

    # 3. Start searching from the year after birth.
    current_year = birth_year + 1

    while True:
        # Define the solar date we are checking for the current year.
        current_solar_date = SolarDate(current_year, birth_month, birth_day)

        # Convert this solar date to its lunar equivalent.
        current_lunar_date = current_solar_date.to_lunar()

        # 4. Check if the lunar month and day match the target lunar birthday.
        if (current_lunar_date.month == target_lunar_month and 
            current_lunar_date.day == target_lunar_day):
            
            # Match found. Print the details as requested.
            print(f"The person's solar birth date is: {birth_year}-{birth_month}-{birth_day}")
            print(f"This corresponds to the lunar month-day: {target_lunar_month}-{target_lunar_day}")
            print(f"The next year when the solar date {birth_month}-{birth_day} aligns with this lunar month-day is:")
            print(current_year)
            break
        
        # If no match, move to the next year.
        current_year += 1
        
        # A safeguard to prevent an infinite loop, although a match is guaranteed.
        if current_year > birth_year + 200:
            print("No matching year found within a 200-year span.")
            break

# Execute the function to find and print the answer.
find_next_matching_birthday()
