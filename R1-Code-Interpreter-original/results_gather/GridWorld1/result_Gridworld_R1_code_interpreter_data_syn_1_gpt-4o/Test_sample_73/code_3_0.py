import json

# Manually constructed path
path = [
    "C5,5", "C5,4", "C4,4", "C3,4", "C3,3", "C2,3", "C2,2", 
    "C2,1", "C3,1", "C4,1", "C4,2", "C4,3", "C4,4", "C4,5", 
    "C3,5", "C2,5", "C1,5"
]

# Output the path as a JSON list
print(json.dumps(path))