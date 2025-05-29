# Convert time to seconds
def time_to_seconds(hours, minutes, seconds):
    return hours * 3600 + minutes * 60 + seconds

# Convert seconds to HH:MM:SS
def seconds_to_time(seconds):
    hours = seconds // 3600
    seconds %= 3600
    minutes = seconds // 60
    seconds %= 60
    return hours, minutes, seconds

# Times in HH:MM:SS
time1 = (23, 12, 56)
time2 = (23, 49, 51)

# Convert both times to seconds
seconds1 = time_to_seconds(*time1)
seconds2 = time_to_seconds(*time2)

# Calculate the difference in seconds
difference_in_seconds = seconds2 - seconds1

# Convert the difference back to HH:MM:SS
difference_time = seconds_to_time(difference_in_seconds)

# Print the result
print(f"{difference_time[0]:02}:{difference_time[1]:02}:{difference_time[2]:02}")