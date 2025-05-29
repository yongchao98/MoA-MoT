import json

# Manually constructed path
path = [
    "C5,5", "C5,4", "C4,4", "C3,4", "C3,3",  # Visit C3,3
    "C3,2", "C2,2",  # Visit C2,2
    "C2,1", "C3,1", "C4,1",  # Visit C4,1
    "C4,2", "C4,3", "C4,4", "C4,5",  # Move to C4,5
    "C5,5", "C5,4",  # Return to C5,4
    "C5,3", "C5,2", "C5,1", "C4,1", "C3,1", "C2,1", "C1,1", "C1,2", "C1,3", "C1,4", "C1,5"  # Visit C1,5
]

# Output the path as a JSON list
print(json.dumps(path))