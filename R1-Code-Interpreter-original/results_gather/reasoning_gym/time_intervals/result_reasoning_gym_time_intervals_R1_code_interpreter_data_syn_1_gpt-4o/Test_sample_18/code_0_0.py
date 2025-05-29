from datetime import datetime

# Define the start and end dates
start_date = datetime(2282, 1, 29)
end_date = datetime(2282, 3, 16)

# Calculate the difference in days
duration = (end_date - start_date).days

# Print the result
print(duration)