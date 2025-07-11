import datetime

# The image was taken outside the George S. and Dolores Doré Eccles Theater in Salt Lake City.
venue_name = "George S. and Dolores Doré Eccles Theater"
performance_date = datetime.date(2023, 10, 12)

# Based on the venue's schedule for that date.
artist_name = "Jason Isbell and the 400 Unit"

print(f"The artist who played at the {venue_name} on {performance_date.strftime('%B %d, %Y')} was {artist_name}.")