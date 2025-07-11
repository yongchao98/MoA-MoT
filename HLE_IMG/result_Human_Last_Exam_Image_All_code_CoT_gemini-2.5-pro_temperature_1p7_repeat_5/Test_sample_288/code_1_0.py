import calendar

# The venue was identified as the Martin Woldson Theater at the Fox in Spokane, WA.
# The task is to find the performer on a specific date.

# Date information
year = 2023
month = 10
day = 12

# Convert month number to month name
month_name = calendar.month_name[month]

# The performer on that date at the identified venue.
# This information was found by searching the venue's past events schedule.
artist_name = "Mat Kearney"

# Print the final answer clearly
print(f"The musical artist that played at this venue on {month_name} {day}, {year} was:")
print(artist_name)