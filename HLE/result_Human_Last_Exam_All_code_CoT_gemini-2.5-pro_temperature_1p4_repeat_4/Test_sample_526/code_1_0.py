import chinesecalendar
from datetime import date

def find_next_matching_birthday():
    """
    Finds the next year after a given birth year where the Solar and Lunar
    birthdays match for a specific date (MM-DD).
    """
    birth_year = 1980
    target_month = 10
    target_day = 1
    
    # We start checking from the year after the birth year.
    # We'll check up to 150 years, which should be more than enough.
    for year in range(birth_year + 1, birth_year + 151):
        try:
            # Create a solar date object for the target day in the current year
            solar_date = date(year, target_month, target_day)
            
            # Convert the solar date to its lunar equivalent
            lunar_date = chinesecalendar.get_lunar_date(solar_date)
            
            # Check if the lunar month and day match the target.
            # We also ensure it's not a leap month for a true match.
            if (lunar_date.month == target_month and 
                lunar_date.day == target_day and 
                not lunar_date.is_leap):
                
                print(f"A person born on 1980-10-01 will have their Solar and Lunar birthdays on the same day (Oct 1st) in the year {year}.")
                print("\nThis means:")
                print(f"Solar Date: {year} - {target_month} - {target_day}")
                print("is the same as")
                # The lunar date includes its own year, which might differ from the solar year, but we'll show all numbers for clarity.
                print(f"Lunar Date: {lunar_date.year} - {lunar_date.month} - {lunar_date.day}")

                return year
        except ValueError:
            # This handles cases where the date might be invalid, although not expected for Oct 1st.
            continue
            
    print("No matching year found within the next 150 years.")
    return None

# Run the function to find the year.
found_year = find_next_matching_birthday()

# The final answer in the requested format is provided separately below.
# Example: <<<2055>>>
if found_year:
    pass # The function already printed the result.
