def find_artist_at_venue():
    """
    This function identifies the musical artist who played at a specific venue on a given date.
    The information was found by first identifying the venue from the image (Martin Woldson Theater at The Fox, Spokane)
    and then looking up its event schedule for October 12, 2023.
    """
    venue_name = "Martin Woldson Theater at The Fox, Spokane"
    performance_date = "October 12, 2023"
    artist_name = "Goo Goo Dolls"

    print(f"The musical artist that played at the {venue_name} on {performance_date} was:")
    print(artist_name)

find_artist_at_venue()