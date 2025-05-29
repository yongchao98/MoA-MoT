from datetime import datetime

# Define the start and end dates
start_date = datetime.strptime("Mon Oct 07 2937", "%a %b %d %Y")
end_date = datetime.strptime("Thu Nov 21 2937", "%a %b %d %Y")

# Calculate the difference in days
difference = end_date - start_date

# Output the number of days
print(difference.days)