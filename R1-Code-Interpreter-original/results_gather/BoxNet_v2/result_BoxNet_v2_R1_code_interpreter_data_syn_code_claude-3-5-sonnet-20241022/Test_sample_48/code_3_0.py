import json

initial_state = {
    "box1": "C5,5", "box2": "C1,3", "box3": "C3,3",
    "box4": "C4,5", "box5": "C4,4"
}

goal_state = {
    "box1": "C2,5", "box2": "C3,4", "box3": "C1,5",
    "box4": "C3,1", "box5": "C1,3"
}

adjacency = {
    "C1,1": ["C1,2", "C2,1"], "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"], "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C2,5"], "C2,1": ["C2,2", "C1,1", "C3,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"], "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"],
    "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"], "C2,5": ["C2,4", "C1,5", "C3,5"],
    "C3,1": ["C3,2", "C2,1", "C4,1"], "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"], "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
    "C3,5": ["C3,4", "C2,5", "C4,5"], "C4,1": ["C4,2", "C3,1", "C5,1"],
    "C4,2": ["C4,1", "C4,3", "C3,2", "C5,2"], "C4,3": ["C4,2", "C4,4", "C3,3", "C5,3"],
    "C4,4": ["C4,3", "C4,5", "C3,4", "C5,4"], "C4,5": ["C4,4", "C3,5", "C5,5"],
    "C5,1": ["C5,2", "C4,1"], "C5,2": ["C5,1", "C5,3", "C4,2"],
    "C5,3": ["C5,2", "C5,4", "C4,3"], "C5,4": ["C5,3", "C5,5", "C4,4"],
    "C5,5": ["C5,4", "C4,5"]
}

# Generate a direct solution
solution = []
current_state = dict(initial_state)
solution.append(dict(current_state))

# Move box1 from C5,5 to C2,5
moves_box1 = [
    {"box1": "C4,5", "box2": "C1,3", "box3": "C3,3", "box4": "C4,5", "box5": "C4,4"},
    {"box1": "C3,5", "box2": "C1,3", "box3": "C3,3", "box4": "C4,5", "box5": "C4,4"},
    {"box1": "C2,5", "box2": "C1,3", "box3": "C3,3", "box4": "C4,5", "box5": "C4,4"}
]

# Move box5 to C1,3 and box2 to C3,4
moves_box5_2 = [
    {"box1": "C2,5", "box2": "C1,3", "box3": "C3,3", "box4": "C4,5", "box5": "C3,4"},
    {"box1": "C2,5", "box2": "C2,3", "box3": "C3,3", "box4": "C4,5", "box5": "C1,3"},
    {"box1": "C2,5", "box2": "C3,4", "box3": "C3,3", "box4": "C4,5", "box5": "C1,3"}
]

# Move box4 to C3,1
moves_box4 = [
    {"box1": "C2,5", "box2": "C3,4", "box3": "C3,3", "box4": "C3,5", "box5": "C1,3"},
    {"box1": "C2,5", "box2": "C3,4", "box3": "C3,3", "box4": "C3,4", "box5": "C1,3"},
    {"box1": "C2,5", "box2": "C3,4", "box3": "C3,3", "box4": "C3,3", "box5": "C1,3"},
    {"box1": "C2,5", "box2": "C3,4", "box3": "C3,3", "box4": "C3,2", "box5": "C1,3"},
    {"box1": "C2,5", "box2": "C3,4", "box3": "C3,3", "box4": "C3,1", "box5": "C1,3"}
]

# Move box3 to C1,5
moves_box3 = [
    {"box1": "C2,5", "box2": "C3,4", "box3": "C2,3", "box4": "C3,1", "box5": "C1,3"},
    {"box1": "C2,5", "box2": "C3,4", "box3": "C1,3", "box4": "C3,1", "box5": "C1,4"},
    {"box1": "C2,5", "box2": "C3,4", "box3": "C1,4", "box4": "C3,1", "box5": "C1,3"},
    {"box1": "C2,5", "box2": "C3,4", "box3": "C1,5", "box4": "C3,1", "box5": "C1,3"}
]

# Combine all moves
solution.extend(moves_box1)
solution.extend(moves_box5_2)
solution.extend(moves_box4)
solution.extend(moves_box3)

print(json.dumps(solution))