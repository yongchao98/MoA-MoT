from datetime import datetime

# Define the two dates
date1_str = "Mon Sep 30 2526"
date2_str = "Wed Nov 27 2526"

# Define the date format
date_format = "%a %b %d %Y"

# Parse the dates into datetime objects
date1 = datetime.strptime(date1_str, date_format)
date2 = datetime.strptime(date2_str, date_format)

# Calculate the difference in days
difference = date2 - date1

# Output the number of days
print(difference.days)