from datetime import datetime

# Define the start and end dates
start_date = datetime.strptime("Fri Feb 20 2652", "%a %b %d %Y")
end_date = datetime.strptime("Tue Feb 24 2652", "%a %b %d %Y")

# Calculate the difference in days
duration = (end_date - start_date).days

# Print the result
print(f"{duration} days")