import json
import sys

def create_solution():
    # Initial state
    states = []
    current = {
        "box1": "C5,5", "box2": "C1,3", "box3": "C3,3", 
        "box4": "C4,5", "box5": "C4,4"
    }
    states.append(dict(current))
    
    # Move box2 to goal
    current["box2"] = "C1,4"; states.append(dict(current))
    current["box2"] = "C2,4"; states.append(dict(current))
    current["box2"] = "C3,4"; states.append(dict(current))
    
    # Move box3 to goal
    current["box3"] = "C2,3"; states.append(dict(current))
    current["box3"] = "C2,4"; states.append(dict(current))
    current["box3"] = "C2,5"; states.append(dict(current))
    current["box3"] = "C1,5"; states.append(dict(current))
    
    # Move box5 to goal
    current["box5"] = "C4,3"; states.append(dict(current))
    current["box5"] = "C3,3"; states.append(dict(current))
    current["box5"] = "C2,3"; states.append(dict(current))
    current["box5"] = "C1,3"; states.append(dict(current))
    
    # Move box4 to goal
    current["box4"] = "C4,4"; states.append(dict(current))
    current["box4"] = "C4,3"; states.append(dict(current))
    current["box4"] = "C4,2"; states.append(dict(current))
    current["box4"] = "C4,1"; states.append(dict(current))
    current["box4"] = "C3,1"; states.append(dict(current))
    
    # Move box1 to goal
    current["box1"] = "C4,5"; states.append(dict(current))
    current["box1"] = "C3,5"; states.append(dict(current))
    current["box1"] = "C2,5"; states.append(dict(current))
    
    return states

# Create solution
solution = create_solution()

# Write solution directly as a string
solution_str = str([{
    "box1": state["box1"],
    "box2": state["box2"],
    "box3": state["box3"],
    "box4": state["box4"],
    "box5": state["box5"]
} for state in solution]).replace("'", '"')

# Print solution directly
sys.stdout.write(solution_str)
sys.stdout.flush()