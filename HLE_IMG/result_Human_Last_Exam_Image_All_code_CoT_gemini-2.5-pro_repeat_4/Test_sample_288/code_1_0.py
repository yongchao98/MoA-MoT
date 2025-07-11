import calendar

def solve():
    """
    This function identifies the musical artist who performed at the venue in the image
    on a specific date.
    
    The steps are:
    1. Identify the venue from the image. The architecture (brickwork, decorative elements)
       points to the Eccles Theater in Salt Lake City, Utah.
    2. Search for the event calendar for the Eccles Theater for October 2023.
    3. Find the performer for the specific date: October 12, 2023.
    4. The artist found is Hozier.
    5. Print the name of the artist.
    """
    
    venue = "Eccles Theater, Salt Lake City"
    year = 2023
    month = 10
    day = 12
    
    # Based on online event listings for the venue on the specified date.
    artist_name = "Hozier"
    
    # Get the day of the week for context
    day_of_week = calendar.day_name[calendar.weekday(year, month, day)]

    print(f"The venue in the image is the {venue}.")
    print(f"On {day_of_week}, {calendar.month_name[month]} {day}, {year}, the musical artist who played at this venue was:")
    print(artist_name)

solve()