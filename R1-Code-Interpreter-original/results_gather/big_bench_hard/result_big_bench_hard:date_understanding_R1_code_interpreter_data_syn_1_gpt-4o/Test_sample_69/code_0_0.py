from datetime import datetime, timedelta

# Jane's birth date
birth_date = datetime(2001, 2, 28)

# Calculate her 16th birthday
sixteenth_birthday = birth_date.replace(year=birth_date.year + 16)

# Print the date in MM/DD/YYYY format
print(sixteenth_birthday.strftime("%m/%d/%Y"))