# Define the start and end times
start_time = "21:39:22"
end_time = "22:25:43"

# Convert time to seconds since the start of the day
def time_to_seconds(time_str):
    h, m, s = map(int, time_str.split(':'))
    return h * 3600 + m * 60 + s

# Calculate the difference in seconds
start_seconds = time_to_seconds(start_time)
end_seconds = time_to_seconds(end_time)
duration_seconds = end_seconds - start_seconds

# Convert the duration back to HH:MM:SS
hours = duration_seconds // 3600
minutes = (duration_seconds % 3600) // 60
seconds = duration_seconds % 60

# Format the result as HH:MM:SS
duration = f"{hours:02}:{minutes:02}:{seconds:02}"
print(duration)