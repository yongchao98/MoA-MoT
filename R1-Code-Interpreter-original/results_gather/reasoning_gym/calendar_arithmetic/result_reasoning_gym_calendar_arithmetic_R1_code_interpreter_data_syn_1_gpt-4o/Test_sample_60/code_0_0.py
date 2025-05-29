import datetime

# Start from the last day of September 2022
date = datetime.date(2022, 9, 30)

# Check each day going backwards until we find a Friday
while date.weekday() != 4:  # 4 corresponds to Friday
    date -= datetime.timedelta(days=1)

# Print the day of the month for the last Friday
print(date.day)