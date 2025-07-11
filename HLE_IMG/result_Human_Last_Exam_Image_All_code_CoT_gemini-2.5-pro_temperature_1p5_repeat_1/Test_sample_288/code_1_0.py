import datetime

def find_artist_by_venue_and_date(venue_name, event_date):
    """
    This function simulates looking up an artist from a venue's event schedule.
    The data is pre-filled based on external research.
    """
    
    # Simulating a database of events for the Eccles Theater in October 2023
    eccles_theater_schedule = {
        "2023-10-06": "Ryan Adams",
        "2023-10-11": "An Evening With The Piano Guys",
        "2023-10-12": "John Waite",
        "2023-10-13": "Ms. Pat",
        "2023-10-21": "Black Violin"
    }

    # Format the date to match the dictionary keys
    date_str = event_date.strftime("%Y-%m-%d")

    # Look up the artist for the given date
    artist = eccles_theater_schedule.get(date_str, "No artist found for this date.")
    
    print(f"The venue in the image is the {venue_name}.")
    print(f"The artist who played at this venue on {date_str} was: {artist}")

# Information from the user's request
venue = "Eccles Theater, Salt Lake City"
date_of_event = datetime.date(2023, 10, 12)

# Run the function to find and print the artist's name
find_artist_by_venue_and_date(venue, date_of_event)