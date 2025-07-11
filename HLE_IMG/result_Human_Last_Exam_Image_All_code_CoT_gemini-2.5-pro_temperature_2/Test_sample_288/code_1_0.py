def find_the_artist():
    """
    This function identifies the venue and then provides the name of the artist
    who performed on the specified date.
    """
    venue_name = "Martin Woldson Theater at The Fox"
    event_date = "October 12, 2023"
    artist_name = "Keb' Mo'"

    print(f"The venue in the image is the {venue_name}.")
    print(f"The artist who performed there on {event_date} was {artist_name}.")

if __name__ == "__main__":
    find_the_artist()