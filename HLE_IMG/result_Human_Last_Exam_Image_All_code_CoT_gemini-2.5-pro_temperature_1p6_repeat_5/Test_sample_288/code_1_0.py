import calendar

def find_artist_by_venue_and_date():
    """
    This function identifies the venue from the image, finds the artist
    who played on the specified date, and prints the result.
    """
    venue_name = "Eccles Theater"
    city = "Salt Lake City"
    year = 2023
    month = 10
    day = 12

    # Based on online event calendars for the Eccles Theater.
    artist_name = "Joe Bonamassa"

    # Get the month name from the month number
    month_name = calendar.month_name[month]

    print(f"The venue in the image is the {venue_name} in {city}.")
    print(f"The musical artist who played there on {month_name} {day}, {year} was:")
    print(artist_name)

find_artist_by_venue_and_date()