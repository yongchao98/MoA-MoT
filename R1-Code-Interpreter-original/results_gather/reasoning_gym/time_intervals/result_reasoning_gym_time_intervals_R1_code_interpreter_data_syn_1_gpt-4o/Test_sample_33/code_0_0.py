from datetime import datetime

# Define the start and end dates
start_date = datetime(2312, 12, 24)
end_date = datetime(2313, 2, 27)

# Calculate the difference in days
duration = (end_date - start_date).days

# Print the result
print(duration)