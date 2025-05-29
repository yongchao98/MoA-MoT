import json

def create_state(box1, box2, box3, box4, box5):
    return {"box1": box1, "box2": box2, "box3": box3, "box4": box4, "box5": box5}

# Generate solution manually with carefully planned moves
solution = [
    # Initial state
    create_state("C5,5", "C1,3", "C3,3", "C4,5", "C4,4"),
    
    # Move box1 to goal (C2,5)
    create_state("C4,5", "C1,3", "C3,3", "C4,5", "C4,4"),
    create_state("C3,5", "C1,3", "C3,3", "C4,5", "C4,4"),
    create_state("C2,5", "C1,3", "C3,3", "C4,5", "C4,4"),
    
    # Move box2 to goal (C3,4)
    create_state("C2,5", "C1,4", "C3,3", "C4,5", "C4,4"),
    create_state("C2,5", "C2,4", "C3,3", "C4,5", "C4,4"),
    create_state("C2,5", "C3,4", "C3,3", "C4,5", "C4,4"),
    
    # Move box3 to goal (C1,5)
    create_state("C2,5", "C3,4", "C2,3", "C4,5", "C4,4"),
    create_state("C2,5", "C3,4", "C1,3", "C4,5", "C4,4"),
    create_state("C2,5", "C3,4", "C1,4", "C4,5", "C4,4"),
    create_state("C2,5", "C3,4", "C1,5", "C4,5", "C4,4"),
    
    # Move box4 to goal (C3,1)
    create_state("C2,5", "C3,4", "C1,5", "C4,4", "C4,3"),
    create_state("C2,5", "C3,4", "C1,5", "C4,3", "C4,2"),
    create_state("C2,5", "C3,4", "C1,5", "C4,2", "C4,1"),
    create_state("C2,5", "C3,4", "C1,5", "C3,2", "C4,1"),
    create_state("C2,5", "C3,4", "C1,5", "C3,1", "C4,1"),
    
    # Move box5 to goal (C1,3)
    create_state("C2,5", "C3,4", "C1,5", "C3,1", "C3,1"),
    create_state("C2,5", "C3,4", "C1,5", "C3,1", "C2,1"),
    create_state("C2,5", "C3,4", "C1,5", "C3,1", "C1,1"),
    create_state("C2,5", "C3,4", "C1,5", "C3,1", "C1,2"),
    create_state("C2,5", "C3,4", "C1,5", "C3,1", "C1,3")
]

print(json.dumps(solution))