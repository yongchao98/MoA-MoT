# Convert time to seconds
def time_to_seconds(hours, minutes, seconds):
    return hours * 3600 + minutes * 60 + seconds

# Convert seconds to HH:MM:SS
def seconds_to_hhmmss(total_seconds):
    hours = total_seconds // 3600
    minutes = (total_seconds % 3600) // 60
    seconds = total_seconds % 60
    return f"{hours:02}:{minutes:02}:{seconds:02}"

# Times given
start_time = (2, 59, 41)  # 02:59:41
end_time = (19, 43, 31)   # 19:43:31

# Convert both times to seconds
start_seconds = time_to_seconds(*start_time)
end_seconds = time_to_seconds(*end_time)

# Calculate the difference in seconds
duration_seconds = end_seconds - start_seconds

# Convert the difference back to HH:MM:SS
duration_hhmmss = seconds_to_hhmmss(duration_seconds)

# Output the result
print(duration_hhmmss)