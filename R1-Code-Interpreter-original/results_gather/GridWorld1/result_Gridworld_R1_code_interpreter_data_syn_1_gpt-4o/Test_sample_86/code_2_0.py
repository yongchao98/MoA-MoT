import json

# Manually constructed path
path = [
    "C4,4", "C4,3", "C5,3",  # Move to C5,3 (goal)
    "C4,3", "C4,2",          # Move to C4,2 (goal)
    "C4,1", "C3,2", "C3,4",  # Move to C3,5 (goal)
    "C3,5", "C2,5",          # Move to C2,5 (goal)
    "C1,5", "C1,4",          # Move to C1,4 (goal)
    "C1,3", "C1,2"           # Move to C1,2 (goal)
]

# Output the path as a JSON list
print(json.dumps(path))