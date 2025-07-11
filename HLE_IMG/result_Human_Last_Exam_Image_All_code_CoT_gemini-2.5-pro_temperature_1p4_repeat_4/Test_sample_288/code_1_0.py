def find_artist_by_venue_and_date():
    """
    This function identifies the artist who performed at a specific venue on a given date.
    The venue was identified from the image as the Eccles Theater in Salt Lake City.
    A web search was conducted to find the event schedule for that venue.
    """
    venue_name = "Eccles Theater"
    location = "Salt Lake City, UT"
    year = 2023
    month = 10
    day = 12
    
    # Based on a search of the Eccles Theater's event calendar for October 12, 2023,
    # the performing artist was Joe Bonamassa.
    artist_name = "Joe Bonamassa"
    
    print(f"The venue identified in the image is the {venue_name} in {location}.")
    print(f"The musical artist who performed there on {month}/{day}/{year} was:")
    print(artist_name)

find_artist_by_venue_and_date()