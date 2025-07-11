import datetime

# The venue was identified as the Paramount Theatre in Denver, CO.
venue_name = "Paramount Theatre"
venue_city = "Denver, CO"

# The specific date requested.
event_date = datetime.date(2023, 10, 12)

# Based on online event schedules for the venue on the given date, the artist was Hozier.
artist_name = "Hozier"

# Print the answer in a full sentence.
print(f"The musical artist that played at the {venue_name} in {venue_city} on {event_date.strftime('%b %d, %Y')} was {artist_name}.")