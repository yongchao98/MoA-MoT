import sys
try:
    from lunardate import LunarDate
except ImportError:
    print("Error: The 'lunardate' library is required. Please install it by running this command in your terminal:")
    print("pip install lunardate")
    sys.exit(1)

def find_next_matching_birthday():
    """
    Finds the next year after 1980 where Solar Oct 1st is also Lunar Oct 1st.
    """
    birth_year = 1980
    
    # Define the target date in both calendars
    target_month = 10
    target_day = 1

    # Start checking from the year after birth
    # A sufficiently large range is used to ensure the date is found.
    for year in range(birth_year + 1, birth_year + 100):
        # Convert the Solar date (YYYY-10-01) for the current year to its Lunar equivalent
        lunar_equiv_date = LunarDate.fromSolarDate(year, target_month, target_day)

        # Check if the Lunar month and day match the target date
        if lunar_equiv_date.month == target_month and lunar_equiv_date.day == target_day:
            # If a match is found, print the final result and the equation
            print(f"The person was born on 1980-10-01.")
            print(f"The next year where the Solar date {target_month}-{target_day} matches the Lunar date {target_month}-{target_day} is:")
            
            # This line outputs each number in the final equation as requested
            print(f"Equation: Solar({year}, {target_month}, {target_day}) == Lunar({lunar_equiv_date.year}, {lunar_equiv_date.month}, {lunar_equiv_date.day})")
            
            print("\nFinal Answer (Year):")
            print(year)
            return

    print(f"No matching year found in the range up to {birth_year + 100}.")

if __name__ == "__main__":
    find_next_matching_birthday()
