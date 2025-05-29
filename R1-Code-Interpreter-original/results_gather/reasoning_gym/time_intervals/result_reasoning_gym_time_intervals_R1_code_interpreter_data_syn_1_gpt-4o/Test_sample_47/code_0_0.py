# Start time
start_hours = 14
start_minutes = 22
start_seconds = 38

# End time
end_hours = 17
end_minutes = 20
end_seconds = 37

# Adjust for borrowing
if end_seconds < start_seconds:
    end_seconds += 60
    end_minutes -= 1

if end_minutes < start_minutes:
    end_minutes += 60
    end_hours -= 1

# Calculate differences
diff_hours = end_hours - start_hours
diff_minutes = end_minutes - start_minutes
diff_seconds = end_seconds - start_seconds

# Print the result
print(f"{diff_hours:02}:{diff_minutes:02}:{diff_seconds:02}")