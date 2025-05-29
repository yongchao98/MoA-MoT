from datetime import datetime, timedelta

# Correct date according to John
correct_date = datetime.strptime("06/19/2019", "%m/%d/%Y")

# Calculate the date one week later
one_week_later = correct_date + timedelta(days=7)

# Format the date in MM/DD/YYYY
formatted_date = one_week_later.strftime("%m/%d/%Y")
print(formatted_date)