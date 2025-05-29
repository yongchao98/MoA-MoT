from datetime import datetime, timedelta

# Start date when Jane bought the eggs
start_date = datetime(2017, 5, 9)

# Number of days until she runs out of eggs
days_to_run_out = 39

# Calculate the date when she runs out of eggs
run_out_date = start_date + timedelta(days=days_to_run_out)

# Add 24 hours to the run out date
final_date = run_out_date + timedelta(days=1)

# Format the final date in MM/DD/YYYY
final_date_str = final_date.strftime("%m/%d/%Y")

print(final_date_str)