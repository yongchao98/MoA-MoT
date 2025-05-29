from datetime import datetime

# Define the two dates
date1 = datetime(2377, 9, 24)
date2 = datetime(2377, 10, 17)

# Calculate the difference in days
duration = (date2 - date1).days

# Print the result
print(duration)