import lunardate
from datetime import date

def find_next_matching_birthday():
    """
    Finds the next year after a given birth year where the Solar date
    corresponds to the same Lunar month and day as the birth date.
    """
    try:
        # Define the person's birth date in the Solar calendar
        birth_year = 1980
        birth_month = 10
        birth_day = 1

        # Convert the Solar birth date to its Lunar equivalent to find the target lunar birthday
        birth_lunar_date = lunardate.LunarDate.fromSolarDate(birth_year, birth_month, birth_day)
        target_lunar_month = birth_lunar_date.month
        target_lunar_day = birth_lunar_date.day

        # Start searching from the year after birth
        current_year = birth_year + 1

        # Loop indefinitely until a match is found or a limit is reached
        while current_year < birth_year + 150: # Set a reasonable search limit
            # Get the Lunar date for October 1st of the current year
            current_lunar_date = lunardate.LunarDate.fromSolarDate(current_year, birth_month, birth_day)

            # Check if the lunar month and day match the target lunar birthday
            if (current_lunar_date.month == target_lunar_month and 
                current_lunar_date.day == target_lunar_day):
                # If they match, we have found the answer
                print(current_year)
                return

            # If no match, move to the next year
            current_year += 1
        
        print("Could not find a matching year within 150 years.")

    except ImportError:
        print("Error: The 'lunardate' library is required.")
        print("Please install it using your package manager, for example: pip install lunardate")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    find_next_matching_birthday()