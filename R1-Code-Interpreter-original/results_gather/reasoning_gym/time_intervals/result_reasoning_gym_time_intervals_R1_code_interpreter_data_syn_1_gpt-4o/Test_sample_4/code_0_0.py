# Define the timestamps
start_time = "01:01:51.896"
end_time = "11:23:19.976"

# Convert the timestamps to total milliseconds
def time_to_milliseconds(time_str):
    hours, minutes, seconds_milliseconds = time_str.split(':')
    seconds, milliseconds = seconds_milliseconds.split('.')
    total_milliseconds = (int(hours) * 3600 * 1000) + (int(minutes) * 60 * 1000) + (int(seconds) * 1000) + int(milliseconds)
    return total_milliseconds

# Calculate the difference in milliseconds
start_ms = time_to_milliseconds(start_time)
end_ms = time_to_milliseconds(end_time)
difference_ms = end_ms - start_ms

# Convert the difference back to HH:MM:SS.mmm format
hours = difference_ms // (3600 * 1000)
difference_ms %= (3600 * 1000)
minutes = difference_ms // (60 * 1000)
difference_ms %= (60 * 1000)
seconds = difference_ms // 1000
milliseconds = difference_ms % 1000

# Format the result
result = f"{hours:02}:{minutes:02}:{seconds:02}.{milliseconds:03}"
print(result)