# Define the start and end times
start_time = (17, 53)  # (hour, minute)
end_time = (17, 59)    # (hour, minute)

# Calculate the difference in minutes
minutes_difference = end_time[1] - start_time[1]

# Express the result in HH:MM format
hours = 0
minutes = minutes_difference

# Print the result
print(f"{hours:02}:{minutes:02}")