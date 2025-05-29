# Define the time slots when Mark is occupied
occupied_times = [
    (5, 10),  # Yoga studio
    (10, 17), # School
    (17, 20), # Rooftop
    (21, 22)  # Garden
]

# Define the bookstore closing time
bookstore_closing_time = 22

# Check for free time slots
free_time_slots = []
for hour in range(24):
    if all(not (start <= hour < end) for start, end in occupied_times) and hour < bookstore_closing_time:
        free_time_slots.append(hour)

# Determine the possible time slot for the bookstore visit
possible_time_slot = None
# Check if there is a continuous free slot from 20 to 21
if 20 in free_time_slots and 21 in free_time_slots:
    possible_time_slot = "8pm to 9pm"

print(possible_time_slot)