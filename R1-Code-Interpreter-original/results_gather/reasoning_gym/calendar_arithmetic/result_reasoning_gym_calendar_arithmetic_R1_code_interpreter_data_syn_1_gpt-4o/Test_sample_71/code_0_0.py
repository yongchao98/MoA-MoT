from datetime import datetime, timedelta

# Define the start and end dates
start_date = datetime(2022, 8, 22)
end_date = datetime(2022, 11, 22)

# Calculate the total number of days including both start and end dates
total_days = (end_date - start_date).days + 1

# Calculate the number of complete weeks and remaining days
complete_weeks = total_days // 7
remaining_days = total_days % 7

# Calculate the number of business days in complete weeks
business_days = complete_weeks * 5

# Calculate the number of business days in the remaining days
for i in range(remaining_days):
    current_day = start_date + timedelta(days=complete_weeks * 7 + i)
    if current_day.weekday() < 5:  # Monday to Friday are business days
        business_days += 1

print(business_days)