import datetime

# The venue was identified as the Eccles Theater in Salt Lake City.
# The event date is provided in the user's request.
venue = "Eccles Theater"
city = "Salt Lake City"
event_date = datetime.date(2023, 10, 12)

# Based on a search of the venue's schedule for that date, the artist was Gryffin.
artist_name = "Gryffin"

# Format the date for printing
formatted_date = event_date.strftime("%B %d, %Y")

print(f"The musical artist that played at the {venue} in {city} on {formatted_date} was {artist_name}.")