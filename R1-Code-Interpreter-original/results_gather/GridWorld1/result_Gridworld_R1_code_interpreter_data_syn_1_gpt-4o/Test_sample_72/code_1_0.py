import json

# Corrected manually constructed path
path = [
    "C5,5", "C5,4",  # Start and visit C5,4
    "C4,4",  # Visit C4,4
    "C3,4", "C3,3",  # Visit C3,3
    "C3,2", "C2,2",  # Visit C2,2
    "C2,1", "C3,1", "C4,1",  # Visit C4,1
    "C4,5", "C3,5", "C2,5", "C1,5"  # Visit C1,5
]

# Output the path as a JSON list
print(json.dumps(path))