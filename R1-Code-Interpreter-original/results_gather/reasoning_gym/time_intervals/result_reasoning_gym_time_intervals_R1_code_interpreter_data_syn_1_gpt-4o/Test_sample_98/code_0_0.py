from datetime import datetime

# Define the two datetime strings
start_date_str = "19 May 2963, 21:51"
end_date_str = "10 August 2963, 14:33"

# Define the format of the datetime strings
date_format = "%d %B %Y, %H:%M"

# Parse the datetime strings into datetime objects
start_date = datetime.strptime(start_date_str, date_format)
end_date = datetime.strptime(end_date_str, date_format)

# Calculate the difference between the two dates
duration = end_date - start_date

# Extract days, hours, minutes, and seconds
days = duration.days
seconds = duration.seconds
hours = seconds // 3600
minutes = (seconds % 3600) // 60
seconds = seconds % 60

# Print the result in the desired format
print(f"{days} days, {hours:02}:{minutes:02}:{seconds:02}")