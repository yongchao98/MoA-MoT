from datetime import datetime, timedelta

# Start date when Jane bought the eggs
start_date = datetime(2017, 5, 9)

# Number of days Jane eats one egg per day
days_to_finish_eggs = 40

# Calculate the date when she finishes the eggs
end_date = start_date + timedelta(days=days_to_finish_eggs)

# Calculate the date for the next day
next_day = end_date + timedelta(days=1)

# Format the date in MM/DD/YYYY
next_day_formatted = next_day.strftime("%m/%d/%Y")

print(next_day_formatted)