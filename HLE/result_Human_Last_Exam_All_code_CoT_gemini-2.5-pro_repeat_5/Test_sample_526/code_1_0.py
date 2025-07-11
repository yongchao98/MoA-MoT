import sys
try:
    import lunardate
    from lunardate import LunarDate
except ImportError:
    print("The 'lunardate' library is required. Please install it by running: pip install lunardate")
    sys.exit(1)

def find_next_matching_birthday():
    """
    Finds the next year after a given birth year where the Solar date MM-DD
    matches the Lunar date MM-DD.
    """
    # The person is born on 1980-10-01.
    birth_year = 1980
    
    # The target birthday is October 1st.
    target_month = 10
    target_day = 1

    print(f"Searching for the next year after {birth_year} where Solar date {target_month:02d}-{target_day:02d} is also the Lunar date {target_month:02d}-{target_day:02d}.")

    # We start searching from the year after birth.
    # A search range of 150 years is more than enough to find the first match.
    for year in range(birth_year + 1, birth_year + 151):
        # Convert the Solar date of the current year to its Lunar equivalent.
        try:
            lunar_equivalent = LunarDate.fromSolarDate(year, target_month, target_day)
        except ValueError as e:
            # This can happen on rare occasions with date inconsistencies, skip the year.
            print(f"Could not process year {year}: {e}")
            continue

        # Check if the Lunar month and day match our target.
        if lunar_equivalent.month == target_month and lunar_equivalent.day == target_day:
            # We found the matching year.
            print(f"\nMatch found!")
            print(f"In the Solar year {year}, the date {year}-{target_month:02d}-{target_day:02d} corresponds to the Lunar date {lunar_equivalent.year}-{lunar_equivalent.month:02d}-{lunar_equivalent.day:02d}.")
            print("\nThe final answer is:")
            print(year)
            return

    print("No match found within the search range.")

if __name__ == "__main__":
    find_next_matching_birthday()
