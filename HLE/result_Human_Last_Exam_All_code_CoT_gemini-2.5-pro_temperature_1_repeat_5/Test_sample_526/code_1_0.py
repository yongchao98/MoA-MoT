import lunardate
from lunardate import Solar

def find_next_matching_birthday():
    """
    Finds the next year after a given birth year where the Solar date (MM-DD)
    corresponds to the same Lunar date (MM-DD) as the birth date.
    """
    # Define the Solar birth date
    birth_year = 1980
    birth_month = 10
    birth_day = 1

    # Convert the Solar birth date to its Lunar equivalent to find our target
    solar_birth_date = Solar(birth_year, birth_month, birth_day)
    lunar_birth_date = solar_birth_date.to_lunar()
    target_lunar_month = lunar_birth_date.month
    target_lunar_day = lunar_birth_date.day

    # Start searching from the year after birth
    current_year = birth_year + 1

    # Loop indefinitely until a match is found
    while True:
        # Define the Solar date for the current year to check
        current_solar_date = Solar(current_year, birth_month, birth_day)
        
        # Convert this Solar date to its Lunar equivalent
        current_lunar_date = current_solar_date.to_lunar()
        
        # Check if the lunar month and day match the target birthday
        if (current_lunar_date.month == target_lunar_month and
            current_lunar_date.day == target_lunar_day):
            
            # If they match, we have found the year. Print it and exit the loop.
            print(current_year)
            break
        
        # If no match, move to the next year
        current_year += 1

if __name__ == '__main__':
    find_next_matching_birthday()