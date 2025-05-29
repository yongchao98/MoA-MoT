from datetime import datetime

# Define the start and end dates
start_date = datetime(2856, 2, 19)
end_date = datetime(2856, 5, 3)

# Calculate the difference in days
duration = (end_date - start_date).days

# Print the result
print(duration)