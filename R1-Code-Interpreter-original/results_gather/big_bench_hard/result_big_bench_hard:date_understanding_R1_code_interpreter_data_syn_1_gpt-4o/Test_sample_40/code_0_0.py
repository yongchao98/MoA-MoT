from datetime import datetime, timedelta

# Jane's birth date
birth_date = datetime(2001, 2, 28)

# Calculate her 16th birthday
sixteenth_birthday = birth_date.replace(year=birth_date.year + 16)

# Calculate the date of the next day
next_day = sixteenth_birthday + timedelta(days=1)

# Format the date in MM/DD/YYYY
next_day_formatted = next_day.strftime("%m/%d/%Y")

print(next_day_formatted)