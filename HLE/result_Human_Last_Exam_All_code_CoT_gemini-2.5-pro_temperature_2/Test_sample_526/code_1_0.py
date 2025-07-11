# To run this code, you may first need to install the lunardate library:
# pip install lunardate

from lunardate import LunarDate
from datetime import date

def find_next_lunar_birthday_match():
    """
    Finds the next year after a given birth year where the Solar and Lunar
    birthdays fall on the same calendar day.
    """
    # Define the original Solar birthday
    birth_year = 1980
    birth_month = 10
    birth_day = 1

    # Convert the Solar birthday to its Lunar equivalent to find our target
    birth_solar_date = date(birth_year, birth_month, birth_day)
    birth_lunar_obj = LunarDate.fromSolarDate(birth_year, birth_month, birth_day)
    target_lunar_month = birth_lunar_obj.month
    target_lunar_day = birth_lunar_obj.day

    print(f"The person was born on the Solar date: {birth_solar_date}")
    print(f"This corresponds to the Lunar date: Month {target_lunar_month}, Day {target_lunar_day}")
    print("-" * 30)
    print(f"Searching for the next year after {birth_year} where {birth_month:02d}-{birth_day:02d} Solar is also {target_lunar_month:02d}-{target_lunar_day:02d} Lunar...")

    # Start searching from the year after birth
    current_year = birth_year + 1

    # Loop indefinitely until a match is found
    while True:
        # Get the Solar date for the current year being checked
        current_solar_date = date(current_year, birth_month, birth_day)

        # Convert it to its Lunar equivalent
        current_lunar_obj = LunarDate.fromSolarDate(current_year, birth_month, birth_day)

        # Check if the lunar month and day match our target
        if (current_lunar_obj.month == target_lunar_month and
            current_lunar_obj.day == target_lunar_day):
            
            print("\nMatch found!")
            print(f"In the year {current_year}, the Solar date {current_solar_date} corresponds to the Lunar date: Month {current_lunar_obj.month}, Day {current_lunar_obj.day}.")
            
            # Print the final answer
            print(f"\nThe next year is: {current_year}")
            break

        # Move to the next year
        current_year += 1
        
        # Safety break to prevent an accidental infinite loop
        if current_year > birth_year + 100:
            print("No match found within 100 years.")
            break

if __name__ == '__main__':
    find_next_lunar_birthday_match()
