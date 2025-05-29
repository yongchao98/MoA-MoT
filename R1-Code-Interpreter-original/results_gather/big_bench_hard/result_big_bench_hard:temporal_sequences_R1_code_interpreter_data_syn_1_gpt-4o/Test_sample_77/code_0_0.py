# Define the time slots when Jennifer was seen
activities = {
    "bike_shop": (6, 9),
    "water_park": (9, 10),
    "electronic_store": (10, 13),
    "yoga_studio": (15, 17)
}

# Define the basketball court open hours
court_open = (0, 17)  # 12am to 5pm

# Check each option to see if it fits the criteria
options = {
    "A": (9, 10),
    "B": (6, 9),
    "C": (15, 17),
    "D": (13, 15)
}

# Function to check if a time slot is free
def is_free(slot, activities):
    for activity in activities.values():
        if not (slot[1] <= activity[0] or slot[0] >= activity[1]):
            return False
    return True

# Find the correct option
for option, time_slot in options.items():
    if is_free(time_slot, activities) and time_slot[1] <= court_open[1]:
        print(option)