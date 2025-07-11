import lunardate

def find_next_matching_birthday():
    """
    Finds the next year after 1980 where the Solar date October 1st
    is also the Lunar date October 1st.
    """
    target_month = 10
    target_day = 1
    
    # The person was born in 1980, so we start checking from the next year.
    year = 1981

    while True:
        # Convert the Solar date (YYYY-10-01) to its Lunar equivalent.
        lunar_date = lunardate.LunarDate.fromSolarDate(year, target_month, target_day)

        # Check if the Lunar month and day match the target month and day.
        if lunar_date.month == target_month and lunar_date.day == target_day:
            # The final equation is Solar YYYY-MM-DD = Lunar YYYY-MM-DD.
            # The numbers in this equation are the year, month, and day.
            # We print the year as the answer.
            print(year)
            break
        
        # If it doesn't match, advance to the next year.
        year += 1

if __name__ == '__main__':
    find_next_matching_birthday()