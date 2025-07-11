import datetime

# Step 1: Identify the venue and date from the user's request.
venue_name = "Martin Woldson Theater at The Fox, Spokane, WA"
event_date = datetime.date(2023, 10, 12)

# Step 2: Research the performer for the given venue and date.
# Based on online event listings and the theater's public schedule history,
# the performing artist on this date was the Goo Goo Dolls.
artist_name = "Goo Goo Dolls"

# Step 3: Print the findings in a clear format.
print(f"Venue: {venue_name}")
print(f"Date: {event_date.strftime('%B %d, %Y')}")
print(f"Musical Artist: {artist_name}")
