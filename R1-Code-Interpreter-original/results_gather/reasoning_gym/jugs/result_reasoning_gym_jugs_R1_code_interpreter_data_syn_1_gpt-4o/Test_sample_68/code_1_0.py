moves = [
    "fill B",
    "pour B->A",
    "empty A",
    "pour B->A",
    "empty A",
    "pour B->A",
    "fill B",
    "pour B->A",
    "empty A",
    "pour B->A",
    "pour B->C"
]

# Print the moves in JSON format
import json
print(json.dumps(moves))