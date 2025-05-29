import json

initial_state = {
    "box1": "C1,4",
    "box2": "C2,8",
    "box3": "C2,2",
    "box4": "C2,1"
}

goal_state = {
    "box1": "C2,8",
    "box2": "C2,2",
    "box3": "C1,2",
    "box4": "C2,5"
}

adjacency = {
    "C1,1": ["C1,2", "C2,1"],
    "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"],
    "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C1,6", "C2,5"],
    "C1,6": ["C1,5", "C1,7", "C2,6"],
    "C1,7": ["C1,6", "C1,8", "C2,7"],
    "C1,8": ["C1,7", "C2,8"],
    "C2,1": ["C2,2", "C1,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3"],
    "C2,4": ["C2,3", "C2,5", "C1,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5"],
    "C2,6": ["C2,5", "C2,7", "C1,6"],
    "C2,7": ["C2,6", "C2,8", "C1,7"],
    "C2,8": ["C2,7", "C1,8"]
}

def verify_move(state1, state2):
    # Count differences
    differences = []
    for box in state1:
        if state1[box] != state2[box]:
            differences.append(box)
    
    # Check only one box moved
    if len(differences) != 1:
        return False
    
    box = differences[0]
    # Check move is to adjacent cell
    if state2[box] not in adjacency[state1[box]]:
        return False
    
    # Check no two boxes in same cell
    if len(set(state2.values())) != len(state2):
        return False
        
    return True

# Create solution path
solution = [
    {"box1": "C1,4", "box2": "C2,8", "box3": "C2,2", "box4": "C2,1"},
    {"box1": "C1,5", "box2": "C2,8", "box3": "C2,2", "box4": "C2,1"},
    {"box1": "C1,6", "box2": "C2,8", "box3": "C2,2", "box4": "C2,1"},
    {"box1": "C1,7", "box2": "C2,8", "box3": "C2,2", "box4": "C2,1"},
    {"box1": "C1,8", "box2": "C2,8", "box3": "C2,2", "box4": "C2,1"},
    {"box1": "C2,8", "box2": "C2,7", "box3": "C2,2", "box4": "C2,1"},
    {"box1": "C2,8", "box2": "C2,6", "box3": "C2,2", "box4": "C2,1"},
    {"box1": "C2,8", "box2": "C2,5", "box3": "C2,2", "box4": "C2,1"},
    {"box1": "C2,8", "box2": "C2,4", "box3": "C2,2", "box4": "C2,1"},
    {"box1": "C2,8", "box2": "C2,3", "box3": "C2,2", "box4": "C2,1"},
    {"box1": "C2,8", "box2": "C2,2", "box3": "C1,2", "box4": "C2,1"},
    {"box1": "C2,8", "box2": "C2,2", "box3": "C1,2", "box4": "C2,2"},
    {"box1": "C2,8", "box2": "C2,2", "box3": "C1,2", "box4": "C2,3"},
    {"box1": "C2,8", "box2": "C2,2", "box3": "C1,2", "box4": "C2,4"},
    {"box1": "C2,8", "box2": "C2,2", "box3": "C1,2", "box4": "C2,5"}
]

# Verify solution
valid = True
if solution[0] != initial_state or solution[-1] != goal_state:
    valid = False
else:
    for i in range(len(solution)-1):
        if not verify_move(solution[i], solution[i+1]):
            valid = False
            break

if valid:
    print(json.dumps(solution))
else:
    print("Invalid solution")