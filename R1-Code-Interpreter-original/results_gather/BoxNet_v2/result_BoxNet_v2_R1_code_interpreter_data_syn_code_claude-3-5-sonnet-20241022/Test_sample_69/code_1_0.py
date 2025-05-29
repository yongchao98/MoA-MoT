import json

def generate_solution():
    states = []
    
    # Initial state
    current = {
        "box1": "C5,5", "box2": "C1,3", "box3": "C3,3",
        "box4": "C4,5", "box5": "C4,4"
    }
    states.append(dict(current))
    
    # 1. First move box1 to goal
    moves_box1 = ["C4,5", "C3,5", "C2,5"]
    for pos in moves_box1:
        current = dict(current)
        current["box1"] = pos
        states.append(dict(current))
    
    # 2. Move box5 to goal
    moves_box5 = ["C3,4", "C2,4", "C1,4", "C1,3"]
    for pos in moves_box5:
        current = dict(current)
        current["box5"] = pos
        states.append(dict(current))
    
    # 3. Move box2 to goal
    moves_box2 = ["C2,3", "C3,3", "C3,4"]
    for pos in moves_box2:
        current = dict(current)
        current["box2"] = pos
        states.append(dict(current))
    
    # 4. Move box4 to goal
    moves_box4 = ["C4,4", "C4,3", "C4,2", "C4,1", "C3,1"]
    for pos in moves_box4:
        current = dict(current)
        current["box4"] = pos
        states.append(dict(current))
    
    # 5. Finally move box3 to goal
    moves_box3 = ["C2,3", "C1,3", "C1,4", "C1,5"]
    for pos in moves_box3:
        current = dict(current)
        current["box3"] = pos
        states.append(dict(current))
    
    return states

# Generate and print solution
solution = generate_solution()
print(json.dumps(solution))