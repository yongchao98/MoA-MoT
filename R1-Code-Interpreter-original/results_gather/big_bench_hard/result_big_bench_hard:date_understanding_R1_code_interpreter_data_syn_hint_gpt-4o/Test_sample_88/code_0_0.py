from datetime import datetime, timedelta

# Start date when Jane bought the eggs
start_date = datetime(2017, 5, 9)

# Number of days she eats one egg per day
days_to_run_out = 40

# Calculate the date when she runs out of eggs
end_date = start_date + timedelta(days=days_to_run_out)

# The date the day after she runs out of eggs
next_day = end_date + timedelta(days=1)

# Format the date in MM/DD/YYYY
next_day_formatted = next_day.strftime("%m/%d/%Y")

print(next_day_formatted)