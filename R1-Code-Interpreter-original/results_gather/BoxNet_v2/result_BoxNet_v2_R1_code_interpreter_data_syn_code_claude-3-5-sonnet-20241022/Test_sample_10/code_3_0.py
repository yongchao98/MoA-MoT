import json
from copy import deepcopy

def get_direct_path():
    # Initial state
    current_state = {
        "box1": "C1,7", "box2": "C2,8", "box3": "C4,5",
        "box4": "C2,4", "box5": "C4,7", "box6": "C3,6"
    }
    
    # Goal state
    goal_state = {
        "box1": "C3,6", "box2": "C4,1", "box3": "C2,7",
        "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"
    }
    
    # Predefined path sequence
    path = [
        # Move box5 (C4,7 -> C2,3)
        {"box1": "C1,7", "box2": "C2,8", "box3": "C4,5", "box4": "C2,4", "box5": "C4,7", "box6": "C3,6"},
        {"box1": "C1,7", "box2": "C2,8", "box3": "C4,5", "box4": "C2,4", "box5": "C4,6", "box6": "C3,6"},
        {"box1": "C1,7", "box2": "C2,8", "box3": "C4,5", "box4": "C2,4", "box5": "C3,6", "box6": "C2,6"},
        {"box1": "C1,7", "box2": "C2,8", "box3": "C4,5", "box4": "C2,4", "box5": "C3,5", "box6": "C2,6"},
        {"box1": "C1,7", "box2": "C2,8", "box3": "C4,5", "box4": "C2,4", "box5": "C3,4", "box6": "C2,6"},
        {"box1": "C1,7", "box2": "C2,8", "box3": "C4,5", "box4": "C2,4", "box5": "C2,4", "box6": "C2,6"},
        {"box1": "C1,7", "box2": "C2,8", "box3": "C4,5", "box4": "C2,3", "box5": "C2,4", "box6": "C2,6"},
        {"box1": "C1,7", "box2": "C2,8", "box3": "C4,5", "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"},
        
        # Move box3 (C4,5 -> C2,7)
        {"box1": "C1,7", "box2": "C2,8", "box3": "C4,6", "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"},
        {"box1": "C1,7", "box2": "C2,8", "box3": "C4,7", "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"},
        {"box1": "C1,7", "box2": "C2,8", "box3": "C3,7", "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"},
        {"box1": "C1,7", "box2": "C2,8", "box3": "C2,7", "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"},
        
        # Move box2 (C2,8 -> C4,1)
        {"box1": "C1,7", "box2": "C3,8", "box3": "C2,7", "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"},
        {"box1": "C1,7", "box2": "C4,8", "box3": "C2,7", "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"},
        {"box1": "C1,7", "box2": "C4,7", "box3": "C2,7", "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"},
        {"box1": "C1,7", "box2": "C4,6", "box3": "C2,7", "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"},
        {"box1": "C1,7", "box2": "C4,5", "box3": "C2,7", "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"},
        {"box1": "C1,7", "box2": "C4,4", "box3": "C2,7", "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"},
        {"box1": "C1,7", "box2": "C4,3", "box3": "C2,7", "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"},
        {"box1": "C1,7", "box2": "C4,2", "box3": "C2,7", "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"},
        {"box1": "C1,7", "box2": "C4,1", "box3": "C2,7", "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"},
        
        # Move box1 (C1,7 -> C3,6)
        {"box1": "C2,7", "box2": "C4,1", "box3": "C2,7", "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"},
        {"box1": "C3,7", "box2": "C4,1", "box3": "C2,7", "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"},
        {"box1": "C3,6", "box2": "C4,1", "box3": "C2,7", "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"}
    ]
    
    print(json.dumps(path))

get_direct_path()