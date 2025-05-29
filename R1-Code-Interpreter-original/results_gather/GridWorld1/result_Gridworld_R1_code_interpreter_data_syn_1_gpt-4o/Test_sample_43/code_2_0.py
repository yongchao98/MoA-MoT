import json

# Define the grid and parameters
initial_position = "C4,4"
goals = {"C2,4", "C1,5", "C3,3", "C4,5"}
obstacles = {"C3,1", "C4,1", "C2,2", "C1,2", "C3,5"}

# Manually construct a valid path
path = [
    "C4,4",  # Start
    "C3,4",  # Move to C3,4
    "C2,4",  # Goal C2,4
    "C1,4",  # Move to C1,4
    "C1,5",  # Goal C1,5
    "C2,5",  # Move to C2,5
    "C3,4",  # Move back to C3,4
    "C3,3",  # Goal C3,3
    "C4,3",  # Move to C4,3
    "C4,4",  # Move back to C4,4
    "C4,5"   # Goal C4,5
]

# Output the path as a JSON list
print(json.dumps(path))