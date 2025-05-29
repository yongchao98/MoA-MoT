import json

# Initial and goal states
initial_state = {
    "box1": "C5,5", "box2": "C1,3", "box3": "C3,3",
    "box4": "C4,5", "box5": "C4,4"
}

goal_state = {
    "box1": "C2,5", "box2": "C3,4", "box3": "C1,5",
    "box4": "C3,1", "box5": "C1,3"
}

# Define the solution step by step
solution = [
    # Initial state
    {"box1": "C5,5", "box2": "C1,3", "box3": "C3,3", "box4": "C4,5", "box5": "C4,4"},
    
    # Move box5 to its goal (C1,3)
    {"box1": "C5,5", "box2": "C1,3", "box3": "C3,3", "box4": "C4,5", "box5": "C3,4"},
    {"box1": "C5,5", "box2": "C1,3", "box3": "C3,3", "box4": "C4,5", "box5": "C2,4"},
    {"box1": "C5,5", "box2": "C2,3", "box3": "C3,3", "box4": "C4,5", "box5": "C2,4"},
    {"box1": "C5,5", "box2": "C2,3", "box3": "C3,3", "box4": "C4,5", "box5": "C1,4"},
    {"box1": "C5,5", "box2": "C2,3", "box3": "C3,3", "box4": "C4,5", "box5": "C1,3"},
    
    # Move box2 to its goal (C3,4)
    {"box1": "C5,5", "box2": "C2,4", "box3": "C3,3", "box4": "C4,5", "box5": "C1,3"},
    {"box1": "C5,5", "box2": "C3,4", "box3": "C3,3", "box4": "C4,5", "box5": "C1,3"},
    
    # Move box3 to its goal (C1,5)
    {"box1": "C5,5", "box2": "C3,4", "box3": "C3,4", "box4": "C4,5", "box5": "C1,3"},
    {"box1": "C5,5", "box2": "C3,4", "box3": "C2,4", "box4": "C4,5", "box5": "C1,3"},
    {"box1": "C5,5", "box2": "C3,4", "box3": "C1,4", "box4": "C4,5", "box5": "C1,3"},
    {"box1": "C5,5", "box2": "C3,4", "box3": "C1,5", "box4": "C4,5", "box5": "C1,3"},
    
    # Move box4 to its goal (C3,1)
    {"box1": "C5,5", "box2": "C3,4", "box3": "C1,5", "box4": "C3,5", "box5": "C1,3"},
    {"box1": "C5,5", "box2": "C3,4", "box3": "C1,5", "box4": "C3,4", "box5": "C1,3"},
    {"box1": "C5,5", "box2": "C3,4", "box3": "C1,5", "box4": "C3,3", "box5": "C1,3"},
    {"box1": "C5,5", "box2": "C3,4", "box3": "C1,5", "box4": "C3,2", "box5": "C1,3"},
    {"box1": "C5,5", "box2": "C3,4", "box3": "C1,5", "box4": "C3,1", "box5": "C1,3"},
    
    # Move box1 to its goal (C2,5)
    {"box1": "C4,5", "box2": "C3,4", "box3": "C1,5", "box4": "C3,1", "box5": "C1,3"},
    {"box1": "C3,5", "box2": "C3,4", "box3": "C1,5", "box4": "C3,1", "box5": "C1,3"},
    {"box1": "C2,5", "box2": "C3,4", "box3": "C1,5", "box4": "C3,1", "box5": "C1,3"}
]

# Verify that each move is valid (only one box moves to an adjacent cell)
def verify_moves():
    for i in range(len(solution) - 1):
        current = solution[i]
        next_state = solution[i + 1]
        
        # Count how many boxes moved
        moves = 0
        moved_box = None
        moved_from = None
        moved_to = None
        
        for box in current:
            if current[box] != next_state[box]:
                moves += 1
                moved_box = box
                moved_from = current[box]
                moved_to = next_state[box]
        
        # Check if exactly one box moved
        if moves != 1:
            print(f"Invalid move at step {i}: {moves} boxes moved")
            return False
        
        # Check if move was to an adjacent cell
        x1, y1 = map(int, moved_from.replace('C', '').split(','))
        x2, y2 = map(int, moved_to.replace('C', '').split(','))
        if abs(x1 - x2) + abs(y1 - y2) != 1:
            print(f"Invalid move at step {i}: non-adjacent move")
            return False
    
    return True

# Verify and print solution
if (solution[0] == initial_state and 
    solution[-1] == goal_state and 
    verify_moves()):
    print(json.dumps(solution))
else:
    print("Invalid solution")