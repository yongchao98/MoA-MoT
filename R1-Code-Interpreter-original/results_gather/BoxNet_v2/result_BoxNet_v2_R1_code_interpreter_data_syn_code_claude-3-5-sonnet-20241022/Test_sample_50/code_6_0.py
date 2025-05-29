import json

def create_fixed_solution():
    # Initial state
    solution = [
        {"box1": "C1,4", "box2": "C2,8", "box3": "C2,2", "box4": "C2,1"},  # Initial
        
        # Move box1 to goal (C2,8)
        {"box1": "C1,5", "box2": "C2,8", "box3": "C2,2", "box4": "C2,1"},
        {"box1": "C1,6", "box2": "C2,8", "box3": "C2,2", "box4": "C2,1"},
        {"box1": "C1,7", "box2": "C2,8", "box3": "C2,2", "box4": "C2,1"},
        {"box1": "C1,8", "box2": "C2,8", "box3": "C2,2", "box4": "C2,1"},
        {"box1": "C2,8", "box2": "C2,7", "box3": "C2,2", "box4": "C2,1"},
        
        # Move box2 to goal (C2,2)
        {"box1": "C2,8", "box2": "C2,6", "box3": "C2,2", "box4": "C2,1"},
        {"box1": "C2,8", "box2": "C2,5", "box3": "C2,2", "box4": "C2,1"},
        {"box1": "C2,8", "box2": "C2,4", "box3": "C2,2", "box4": "C2,1"},
        {"box1": "C2,8", "box2": "C2,3", "box3": "C2,2", "box4": "C2,1"},
        {"box1": "C2,8", "box2": "C2,2", "box3": "C2,3", "box4": "C2,1"},
        
        # Move box3 to goal (C1,2)
        {"box1": "C2,8", "box2": "C2,2", "box3": "C1,3", "box4": "C2,1"},
        {"box1": "C2,8", "box2": "C2,2", "box3": "C1,2", "box4": "C2,1"},
        
        # Move box4 to goal (C2,5)
        {"box1": "C2,8", "box2": "C2,2", "box3": "C1,2", "box4": "C2,2"},
        {"box1": "C2,8", "box2": "C2,2", "box3": "C1,2", "box4": "C2,3"},
        {"box1": "C2,8", "box2": "C2,2", "box3": "C1,2", "box4": "C2,4"},
        {"box1": "C2,8", "box2": "C2,2", "box3": "C1,2", "box4": "C2,5"}
    ]
    return solution

# Generate and print solution
solution = create_fixed_solution()
print(json.dumps(solution))