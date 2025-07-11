import calendar

# Step 1: Identify the venue.
# Based on the architecture (red brick with herringbone pattern, tan sections, quatrefoil emblem),
# the venue is identified as the George S. and Dolores Doré Eccles Theater in Salt Lake City, Utah.
venue_name = "Eccles Theater"
city = "Salt Lake City"

# Step 2: Identify the date from the user's request.
year = 2023
month = 10
day = 12

# Step 3: Find the performer on that specific date at the venue.
# A search of the Eccles Theater's schedule for October 2023 reveals the musical artist.
artist_name = "The National"

# Step 4: Print the answer in a clear statement.
month_name = calendar.month_name[month]
print(f"The musical artist that played at the {venue_name} in {city} on {month_name} {day}, {year} was {artist_name}.")