# Define the vehicles and their positions
vehicles = ["limousine", "truck", "minivan", "tractor", "hatchback", "bus", "convertible"]

# Define the conditions
conditions = {
    "convertible": 2,
    "hatchback": 3,
    "tractor": 4,
    "minivan": None,
    "truck": None,
    "limousine": None,
    "bus": None
}

# Determine the order based on conditions
order = [None] * 7
order[1] = "convertible"
order[2] = "hatchback"
order[3] = "tractor"

# Place minivan, truck, and limousine based on their relative order
order[0] = "minivan"
order[4] = "truck"
order[5] = "limousine"
order[6] = "bus"

# Check the oldest vehicle
oldest_vehicle = order[0]
print(oldest_vehicle)