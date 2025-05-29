sequence = [9, 2, 2, 1, 0, 0, 0]
counts = {}

for num in sequence:
    if num in counts:
        counts[num] += 1
    else:
        counts[num] = 1

# Find the next number based on the pattern
last_number = sequence[-1]
next_number = None

# Check if the last number has appeared the required number of times
if counts[last_number] < 3:
    next_number = last_number
else:
    # Find the next number that should appear
    for num in range(last_number - 1, -1, -1):
        if num not in counts or counts[num] < num:
            next_number = num
            break

print(next_number)