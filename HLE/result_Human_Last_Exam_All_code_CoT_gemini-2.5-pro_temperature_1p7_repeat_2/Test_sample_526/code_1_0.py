import lunardate
from datetime import date

def find_next_matching_birthday():
    """
    Finds the next year after 1980 where the Solar date 10-01
    is also the Lunar date 10-01.
    """
    # The person's Solar birthday month and day.
    solar_month = 10
    solar_day = 1

    # We are looking for a year where the Lunar birthday is the same MM-DD.
    target_lunar_month = 10
    target_lunar_day = 1

    # Start checking from the year after the birth year.
    start_year = 1980 + 1
    # Set a reasonable end year to avoid an infinite loop.
    end_year = 2100

    found_year = None
    for year in range(start_year, end_year + 1):
        try:
            # Convert the Solar date (YYYY-10-01) of the current iteration year to Lunar.
            lunar_date_obj = lunardate.LunarDate.fromSolarDate(year, solar_month, solar_day)

            # Check if the resulting lunar month and day match the target.
            if lunar_date_obj.month == target_lunar_month and lunar_date_obj.day == target_lunar_day:
                found_year = year
                break  # Exit the loop as we found the first matching year.
        except (ValueError, TypeError) as e:
            # Handle potential errors with date conversion if any.
            print(f"Could not process year {year}: {e}")
            continue

    if found_year:
        print(f"The next year since their birth that their Lunar calendar birthday (10-01) is the same as their Solar Calendar's birthday (10-01) is: {found_year}")
        print(f"<<<{found_year}>>>")
    else:
        print(f"No matching year was found before {end_year}.")
        print("<<<Not Found>>>")

find_next_matching_birthday()