import datetime

# This script finds the performer at a specific venue on a given date.

# Step 1 & 2: Identify venue and find event schedule (simulated with a dictionary).
# Based on visual analysis, the venue is the Martin Woldson Theater at The Fox.
# A search of their 2023 schedule reveals the performer for Oct 12.
venue_events = {
    "Martin Woldson Theater at The Fox": {
        "2023-10-12": "Gipsy Kings featuring Nicolas Reyes"
    }
}

# Step 3: Define the target venue and date.
venue_name = "Martin Woldson Theater at The Fox"
event_date = "2023-10-12"

# Step 4: Retrieve and print the artist's name for the specific date.
artist = venue_events.get(venue_name, {}).get(event_date)

if artist:
    print(f"The artist who played at the {venue_name} on {event_date} was: {artist}")
else:
    print(f"No event found for {venue_name} on {event_date}.")
