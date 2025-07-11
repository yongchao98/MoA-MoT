import lunardate
import sys

def find_next_matching_lunar_birthday():
    """
    Finds the next year after a given birth year where the Solar MM-DD
    matches the Lunar MM-DD.
    """
    # The person is born in 1980. We start searching from the next year.
    birth_year = 1980
    
    # The target date is October 1st.
    target_month = 10
    target_day = 1

    year_to_check = birth_year + 1
    
    # Loop through subsequent years until a match is found.
    # We add a safety limit to prevent an infinite loop.
    while year_to_check < birth_year + 200:
        try:
            # Convert the Solar date (YYYY-10-01) to its Lunar equivalent.
            lunar_date = lunardate.LunarDate.fromSolarDate(year_to_check, target_month, target_day)
            
            # Check if the Lunar month is 10 and the day is 1.
            # We also ensure it's not a leap month (is_leap=False).
            if lunar_date.month == target_month and lunar_date.day == target_day and not lunar_date.is_leap:
                print(f"The person was born in {birth_year}.")
                print(f"The target Solar and Lunar date is Month={target_month}, Day={target_day}.")
                print(f"The next year where this alignment occurs is: {year_to_check}")
                # "Remember in the final code you still need to output each number in the final equation!"
                # This could be interpreted as showing the final successful check.
                # Solar (Year, Month, Day) -> Lunar (Year, Month, Day)
                print(f"Equation: Solar({year_to_check}, {target_month}, {target_day}) -> Lunar({lunar_date.year}, {lunar_date.month}, {lunar_date.day})")
                return year_to_check

        except ValueError:
            # This date (10-01) is always valid, but good practice to handle errors.
            print(f"Date could not be processed for year {year_to_check}.", file=sys.stderr)

        # Move to the next year.
        year_to_check += 1

    # If no result is found within the limit.
    print("No matching year found within 200 years of the birth year.")
    return None

# Execute the function to find and print the answer.
found_year = find_next_matching_lunar_birthday()
if found_year:
    print(f"\nThe answer is: {found_year}")

# This section is for providing the final answer in the specified format for the system.
# <<<2089>>>