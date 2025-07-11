import datetime

# Time the bottles were found open
time_found_open_str = "1:00 PM"

# Time the bottles were left open from
time_left_open_str = "7:00 AM"

# Parse the times
time_found_open = datetime.datetime.strptime(time_found_open_str, "%I:%M %p").time()
time_left_open = datetime.datetime.strptime(time_left_open_str, "%I:%M %p").time()

# Create datetime objects for calculation (date is arbitrary, we only need the time difference)
today = datetime.date.today()
datetime1 = datetime.datetime.combine(today, time_left_open)
datetime2 = datetime.datetime.combine(today, time_found_open)

# Calculate the difference
exposure_duration = datetime2 - datetime1
exposure_hours = exposure_duration.total_seconds() / 3600

# Print key experimental parameters
print(f"Analysis of Experimental Parameters:")
print(f"Autoclave Temperature: 121 degrees Celsius")
print(f"Autoclave Time: 25 minutes")
print(f"Media Exposure Duration to Air: The PDA bottles were left open for {int(exposure_hours)} hours.")
print(f"Incubation Temperature: 25.5 degrees Celsius")
print(f"Incubation Duration: 5 days")
