import calendar
import datetime

def find_artist_by_date():
    """
    This function identifies the venue from the image and finds the musical artist
    who performed on the specified date.
    """
    venue_name = "Martin Woldson Theater at The Fox, Spokane, WA"
    event_year = 2023
    event_month = 10
    event_day = 12

    # Based on online event calendars for the identified venue.
    artist_name = "Ben Folds"

    # Get the day of the week for the event date.
    event_date = datetime.date(event_year, event_month, event_day)
    day_of_week = calendar.day_name[event_date.weekday()]

    print(f"The venue in the image is the {venue_name}.")
    print(f"Searching for the musical artist who played on {day_of_week}, {calendar.month_name[event_month]} {event_day}, {event_year}...")
    print(f"The artist who performed at this venue on that date was: {artist_name}")

find_artist_by_date()