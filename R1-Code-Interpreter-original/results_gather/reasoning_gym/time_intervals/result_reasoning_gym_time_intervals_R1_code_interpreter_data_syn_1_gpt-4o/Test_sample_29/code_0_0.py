from datetime import datetime

# Define the start and end dates
start_date = datetime(2699, 2, 7)
end_date = datetime(2699, 3, 6)

# Calculate the difference in days
difference = (end_date - start_date).days

# Print the result
print(difference)