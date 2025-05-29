# Define the start and end times
start_time = (15, 32, 18)  # 15:32:18
end_time = (23, 2, 38)    # 23:02:38

# Convert both times to seconds since the start of the day
start_seconds = start_time[0] * 3600 + start_time[1] * 60 + start_time[2]
end_seconds = end_time[0] * 3600 + end_time[1] * 60 + end_time[2]

# Calculate the difference in seconds
elapsed_seconds = end_seconds - start_seconds

# Convert the difference back to hours, minutes, and seconds
hours = elapsed_seconds // 3600
minutes = (elapsed_seconds % 3600) // 60
seconds = elapsed_seconds % 60

# Print the result in HH:MM:SS format
print(f"{hours:02}:{minutes:02}:{seconds:02}")