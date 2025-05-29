# Initial and goal states
initial_state = {
    "box1": "C1,3",
    "box2": "C2,4",
    "box3": "C1,6",
    "box4": "C2,2",
    "box5": "C1,4"
}

goal_state = {
    "box1": "C2,1",
    "box2": "C1,2",
    "box3": "C1,5",
    "box4": "C2,2",
    "box5": "C2,6"
}

# Adjacency list
adjacency = {
    "C1,1": ["C1,2", "C2,1"],
    "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"],
    "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C1,6", "C2,5"],
    "C1,6": ["C1,5", "C2,6"],
    "C2,1": ["C2,2", "C1,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3"],
    "C2,4": ["C2,3", "C2,5", "C1,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5"],
    "C2,6": ["C2,5", "C1,6"]
}

# Plan generation
plan = [initial_state.copy()]

# Move box1 from C1,3 to C2,1
plan.append({"box1": "C1,2", "box2": "C2,4", "box3": "C1,6", "box4": "C2,2", "box5": "C1,4"})
plan.append({"box1": "C1,1", "box2": "C2,4", "box3": "C1,6", "box4": "C2,2", "box5": "C1,4"})
plan.append({"box1": "C2,1", "box2": "C2,4", "box3": "C1,6", "box4": "C2,2", "box5": "C1,4"})

# Move box2 from C2,4 to C1,2
plan.append({"box1": "C2,1", "box2": "C2,3", "box3": "C1,6", "box4": "C2,2", "box5": "C1,4"})
plan.append({"box1": "C2,1", "box2": "C2,2", "box3": "C1,6", "box4": "C2,2", "box5": "C1,4"})
plan.append({"box1": "C2,1", "box2": "C1,2", "box3": "C1,6", "box4": "C2,2", "box5": "C1,4"})

# Move box3 from C1,6 to C1,5
plan.append({"box1": "C2,1", "box2": "C1,2", "box3": "C1,5", "box4": "C2,2", "box5": "C1,4"})

# Move box5 from C1,4 to C2,6
plan.append({"box1": "C2,1", "box2": "C1,2", "box3": "C1,5", "box4": "C2,2", "box5": "C1,5"})
plan.append({"box1": "C2,1", "box2": "C1,2", "box3": "C1,5", "box4": "C2,2", "box5": "C1,6"})
plan.append({"box1": "C2,1", "box2": "C1,2", "box3": "C1,5", "box4": "C2,2", "box5": "C2,6"})

# Print the plan
print(plan)