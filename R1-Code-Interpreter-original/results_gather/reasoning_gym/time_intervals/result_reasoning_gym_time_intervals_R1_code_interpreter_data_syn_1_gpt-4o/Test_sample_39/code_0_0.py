# Define the start and end times
start_time = "16:39:51.271"
end_time = "16:44:29.429"

# Convert time to hours, minutes, seconds, and milliseconds
def time_to_milliseconds(time_str):
    hours, minutes, seconds_milliseconds = time_str.split(':')
    seconds, milliseconds = seconds_milliseconds.split('.')
    total_milliseconds = (int(hours) * 3600 * 1000) + (int(minutes) * 60 * 1000) + (int(seconds) * 1000) + int(milliseconds)
    return total_milliseconds

# Calculate the difference in milliseconds
start_milliseconds = time_to_milliseconds(start_time)
end_milliseconds = time_to_milliseconds(end_time)
duration_milliseconds = end_milliseconds - start_milliseconds

# Convert milliseconds back to HH:MM:SS.mmm
hours = duration_milliseconds // (3600 * 1000)
duration_milliseconds %= (3600 * 1000)
minutes = duration_milliseconds // (60 * 1000)
duration_milliseconds %= (60 * 1000)
seconds = duration_milliseconds // 1000
milliseconds = duration_milliseconds % 1000

# Format the result
duration = f"{hours:02}:{minutes:02}:{seconds:02}.{milliseconds:03}"
print(duration)