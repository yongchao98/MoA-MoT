import sys

def solve():
    """
    Finds the next year after 1980 where the Solar and Lunar birthdays match for a person born on 1980-10-01.
    """
    try:
        from sxtwl import fromSolar, fromLunar
    except ImportError:
        print("Error: The 'sxtwl' library is required. Please install it by running: pip install sxtwl", file=sys.stderr)
        sys.exit(1)

    # The person's Solar birth date
    birth_year = 1980
    birth_month = 10
    birth_day = 1

    # Step 1: Convert the Solar birth date to its Lunar equivalent.
    solar_birth_obj = fromSolar(birth_year, birth_month, birth_day)
    lunar_birth_month = solar_birth_obj.month
    lunar_birth_day = solar_birth_obj.day
    # The Lunar date is not in a leap month, which simplifies the search.
    is_leap = solar_birth_obj.isLeap

    print(f"Solar Birthday: {birth_year}-{birth_month:02d}-{birth_day:02d}")
    print(f"Equivalent Lunar Birthday: Month {lunar_birth_month}, Day {lunar_birth_day}")
    print("-" * 30)

    # Step 2: Loop through years starting from the year after birth.
    search_year = birth_year + 1
    # Set a reasonable upper limit for the search.
    end_year = birth_year + 100 

    while search_year <= end_year:
        try:
            # Step 3: Convert the Lunar birthday for the current search year back to Solar.
            solar_equivalent_obj = fromLunar(search_year, lunar_birth_month, lunar_birth_day, is_leap)
            
            # Step 4: Check if the resulting Solar month and day match the original.
            if solar_equivalent_obj.month == birth_month and solar_equivalent_obj.day == birth_day:
                
                # Step 5: Match found. Print the details and the final answer.
                print(f"Searching in year: {search_year}")
                print(f"The Lunar date {lunar_birth_month:02d}-{lunar_birth_day:02d} corresponds to the Solar date: {solar_equivalent_obj.year}-{solar_equivalent_obj.month:02d}-{solar_equivalent_obj.day:02d}")
                print(f"This matches the original birthday's month and day ({birth_month:02d}-{birth_day:02d}).")
                print("\nFinal Answer:")
                print(f"{search_year}")
                
                # Return the final answer in the required format
                return search_year

        except Exception:
            # This can happen if a lunar date (e.g., month 8, day 30) doesn't exist in a particular year.
            # We can safely ignore these years and continue searching.
            pass
            
        search_year += 1

    print(f"No match found between {birth_year + 1} and {end_year}.")
    return None

# Execute the function to find and print the answer.
found_year = solve()

if found_year:
    print(f"<<<{found_year}>>>")