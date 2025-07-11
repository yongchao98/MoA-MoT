import calendar

# Step 1: Identify the venue and date from the user's request.
# The image shows the Eccles Theater in Salt Lake City.
venue = "Eccles Theater"
year = 2023
month = 10
day = 12

# Step 2: Research the performer for the given date and venue.
# A search of the venue's schedule for October 12, 2023, shows a concert by John Mayer.
artist = "John Mayer"

# Step 3: Print the result in a clear, formatted sentence.
# Use the calendar module to get the month name.
month_name = calendar.month_name[month]

print(f"The musical artist who played at the {venue} on {month_name} {day}, {year} was:")
print(artist)