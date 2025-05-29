import json

def create_state(box1, box2, box3, box4, box5):
    return {
        "box1": box1,
        "box2": box2,
        "box3": box3,
        "box4": box4,
        "box5": box5
    }

# Create minimal solution with careful planning
solution = [
    # Initial state
    create_state("C5,5", "C1,3", "C3,3", "C4,5", "C4,4"),
    
    # Move box1 to C2,5
    create_state("C4,5", "C1,3", "C3,3", "C4,5", "C4,4"),
    create_state("C3,5", "C1,3", "C3,3", "C4,5", "C4,4"),
    create_state("C2,5", "C1,3", "C3,3", "C4,5", "C4,4"),
    
    # Move box4 to C3,1
    create_state("C2,5", "C1,3", "C3,3", "C4,4", "C4,3"),
    create_state("C2,5", "C1,3", "C3,3", "C4,3", "C4,2"),
    create_state("C2,5", "C1,3", "C3,3", "C4,2", "C4,1"),
    create_state("C2,5", "C1,3", "C3,3", "C4,1", "C3,1"),
    create_state("C2,5", "C1,3", "C3,3", "C3,1", "C3,2"),
    
    # Move box2 to C3,4
    create_state("C2,5", "C2,3", "C3,3", "C3,1", "C3,2"),
    create_state("C2,5", "C3,3", "C3,2", "C3,1", "C3,3"),
    create_state("C2,5", "C3,4", "C3,2", "C3,1", "C3,3"),
    
    # Move box3 to C1,5
    create_state("C2,5", "C3,4", "C2,2", "C3,1", "C3,3"),
    create_state("C2,5", "C3,4", "C1,2", "C3,1", "C3,3"),
    create_state("C2,5", "C3,4", "C1,3", "C3,1", "C3,3"),
    create_state("C2,5", "C3,4", "C1,4", "C3,1", "C3,3"),
    create_state("C2,5", "C3,4", "C1,5", "C3,1", "C3,3"),
    
    # Move box5 to C1,3
    create_state("C2,5", "C3,4", "C1,5", "C3,1", "C2,3"),
    create_state("C2,5", "C3,4", "C1,5", "C3,1", "C1,3")
]

print(json.dumps(solution))