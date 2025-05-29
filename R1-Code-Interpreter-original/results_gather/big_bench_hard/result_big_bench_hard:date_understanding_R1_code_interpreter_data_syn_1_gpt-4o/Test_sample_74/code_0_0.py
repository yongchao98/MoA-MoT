from datetime import datetime, timedelta

# Define the current date
current_date = datetime(2003, 9, 7)

# Calculate the date one year ago
one_year_ago = current_date.replace(year=current_date.year - 1)

# Format the date in MM/DD/YYYY
formatted_date = one_year_ago.strftime("%m/%d/%Y")
print(formatted_date)