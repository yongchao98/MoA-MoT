from datetime import datetime, timedelta

# Jane's birth date
birth_date = datetime(2000, 2, 29)

# Calculate her 16th birthday
sixteenth_birthday = birth_date.replace(year=birth_date.year + 16)

# Calculate the next day
next_day = sixteenth_birthday + timedelta(days=1)

# Format the date as MM/DD/YYYY
formatted_date = next_day.strftime("%m/%d/%Y")

print(formatted_date)