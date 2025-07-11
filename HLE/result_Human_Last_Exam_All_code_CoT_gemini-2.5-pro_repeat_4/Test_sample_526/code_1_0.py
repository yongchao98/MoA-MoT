# The user might need to install the 'lunardate' library first by running:
# pip install lunardate

from lunardate import LunarDate

def find_next_matching_birthday_year():
    """
    Finds the next year after 1980 where the Solar date October 1st
    is also the Lunar date October 1st.
    """
    # Start checking from the year after the person's birth year (1980).
    year = 1981

    # Loop indefinitely until the matching year is found.
    while True:
        try:
            # Convert the Solar date (year, 10, 1) to its Lunar equivalent.
            lunar_date = LunarDate.fromSolarDate(year, 10, 1)

            # Check if the Lunar date is also October 1st.
            if lunar_date.month == 10 and lunar_date.day == 1:
                # If a match is found, print the year and exit the loop.
                print(year)
                break
        except (ValueError, OverflowError):
            # In case of any conversion error for a given year,
            # we simply skip to the next year.
            pass
        
        # Move to the next year.
        year += 1

if __name__ == "__main__":
    find_next_matching_birthday_year()