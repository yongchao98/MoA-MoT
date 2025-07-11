import sxtwl

def find_next_lunar_birthday_match():
    """
    Finds the next year after a given birth year where the solar birthday
    (e.g., Oct 1st) falls on the same lunar date as the original birth date.
    """
    # Define the person's solar birth date
    birth_year = 1980
    birth_month = 10
    birth_day = 1

    # Step 1: Find the corresponding lunar birthday from the birth date.
    # This becomes the target lunar date we need to find in a future year.
    birth_date_obj = sxtwl.fromSolar(birth_year, birth_month, birth_day)
    target_lunar_month = birth_date_obj.getLunarMonth()
    target_lunar_day = birth_date_obj.getLunarDay()
    
    # We must also note if the original birthday was in a leap month (it wasn't).
    # We will only match with non-leap months.
    target_is_leap = birth_date_obj.isLunarLeap()

    # Step 2: Search for the next year where the solar birthday aligns with the lunar birthday.
    # We start from the year after birth.
    search_year = birth_year + 1
    
    # Loop through a reasonable number of years (e.g., 100 years).
    while search_year < birth_year + 100:
        # For the current year being checked, get the lunar date for the solar birthday.
        current_date_obj = sxtwl.fromSolar(search_year, birth_month, birth_day)
        
        current_lunar_month = current_date_obj.getLunarMonth()
        current_lunar_day = current_date_obj.getLunarDay()
        current_is_leap = current_date_obj.isLunarLeap()
        
        # Step 3: Check if the lunar date for this year's birthday matches the original lunar birthday.
        if (current_lunar_month == target_lunar_month and 
            current_lunar_day == target_lunar_day and 
            current_is_leap == target_is_leap):
            
            # If a match is found, print the details and the final answer.
            print(f"A person born on the solar date {birth_year}-{birth_month}-{birth_day} has a lunar birthday of month {target_lunar_month}, day {target_lunar_day}.")
            print(f"We are looking for the next year when the solar date {birth_month}-{birth_day} again falls on this specific lunar birthday.")
            print("\n--- Match Found ---")
            print(f"In the year {search_year}, the solar date {search_year}-{birth_month}-{birth_day} corresponds to the lunar date of month {current_lunar_month}, day {current_lunar_day}.")
            print("\nThe final answer is the year:")
            print(search_year)
            return search_year
            
        search_year += 1
        
    print("No match found within 100 years.")
    return None

# Run the function to find and print the answer.
find_next_lunar_birthday_match()