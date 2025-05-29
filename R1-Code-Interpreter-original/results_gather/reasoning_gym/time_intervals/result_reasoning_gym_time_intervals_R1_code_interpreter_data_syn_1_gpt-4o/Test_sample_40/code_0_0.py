# Start and end times
start_hour, start_minute = 19, 47
end_hour, end_minute = 21, 21

# Calculate the difference
hours_difference = end_hour - start_hour
minutes_difference = end_minute - start_minute

# Adjust if minutes difference is negative
if minutes_difference < 0:
    hours_difference -= 1
    minutes_difference += 60

# Print the result in HH:MM format
print(f"{hours_difference:02}:{minutes_difference:02}")