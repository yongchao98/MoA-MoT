import json

# Create the solution step by step
solution = []

# Initial state
solution.append({"box1": "C1,4", "box2": "C2,8", "box3": "C2,2", "box4": "C2,1"})

# Move box3 to goal (C1,2)
solution.append({"box1": "C1,4", "box2": "C2,8", "box3": "C1,2", "box4": "C2,1"})

# Move box2 to goal (C2,2)
solution.append({"box1": "C1,4", "box2": "C2,7", "box3": "C1,2", "box4": "C2,1"})
solution.append({"box1": "C1,4", "box2": "C2,6", "box3": "C1,2", "box4": "C2,1"})
solution.append({"box1": "C1,4", "box2": "C2,5", "box3": "C1,2", "box4": "C2,1"})
solution.append({"box1": "C1,4", "box2": "C2,4", "box3": "C1,2", "box4": "C2,1"})
solution.append({"box1": "C1,4", "box2": "C2,3", "box3": "C1,2", "box4": "C2,1"})
solution.append({"box1": "C1,4", "box2": "C2,2", "box3": "C1,2", "box4": "C2,1"})

# Move box4 to goal (C2,5)
solution.append({"box1": "C1,4", "box2": "C2,2", "box3": "C1,2", "box4": "C2,2"})
solution.append({"box1": "C1,4", "box2": "C2,2", "box3": "C1,2", "box4": "C2,3"})
solution.append({"box1": "C1,4", "box2": "C2,2", "box3": "C1,2", "box4": "C2,4"})
solution.append({"box1": "C1,4", "box2": "C2,2", "box3": "C1,2", "box4": "C2,5"})

# Move box1 to goal (C2,8)
solution.append({"box1": "C1,5", "box2": "C2,2", "box3": "C1,2", "box4": "C2,5"})
solution.append({"box1": "C1,6", "box2": "C2,2", "box3": "C1,2", "box4": "C2,5"})
solution.append({"box1": "C1,7", "box2": "C2,2", "box3": "C1,2", "box4": "C2,5"})
solution.append({"box1": "C1,8", "box2": "C2,2", "box3": "C1,2", "box4": "C2,5"})
solution.append({"box1": "C2,8", "box2": "C2,2", "box3": "C1,2", "box4": "C2,5"})

# Convert to string with fixed formatting
result = "["
for i, state in enumerate(solution):
    if i > 0:
        result += ", "
    result += json.dumps(state)
result += "]"

# Print with proper wrapping
print("<<<" + result + ">>>")