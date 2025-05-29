from datetime import datetime

# Define the start and end dates
start_date = datetime(2560, 3, 11)
end_date = datetime(2560, 3, 22)

# Calculate the difference in days
delta = end_date - start_date

# Output the number of days
print(f"{delta.days} days")