# Define start and end times in HH:MM:SS
start_time = (12, 32, 23)
end_time = (14, 58, 48)

# Convert start and end times to total seconds
start_seconds = start_time[0] * 3600 + start_time[1] * 60 + start_time[2]
end_seconds = end_time[0] * 3600 + end_time[1] * 60 + end_time[2]

# Calculate the duration in seconds
duration_seconds = end_seconds - start_seconds

# Convert the duration back to HH:MM:SS
hours = duration_seconds // 3600
minutes = (duration_seconds % 3600) // 60
seconds = duration_seconds % 60

# Print the result in HH:MM:SS format
print(f"{hours:02}:{minutes:02}:{seconds:02}")