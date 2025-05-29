from datetime import datetime

# Define the start and end dates
start_date = datetime(2322, 2, 25)
end_date = datetime(2322, 3, 29)

# Calculate the difference in days
duration = (end_date - start_date).days

# Print the result
print(duration)