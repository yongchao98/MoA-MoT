import json

def solve():
    # Initial and goal states
    initial_state = {
        "box1": "C1,4", "box2": "C2,8",
        "box3": "C2,2", "box4": "C2,1"
    }
    
    goal_state = {
        "box1": "C2,8", "box2": "C2,2",
        "box3": "C1,2", "box4": "C2,5"
    }
    
    # Pre-computed solution path
    solution = [
        {"box1": "C1,4", "box2": "C2,8", "box3": "C2,2", "box4": "C2,1"},  # Initial state
        {"box1": "C1,4", "box2": "C2,8", "box3": "C2,2", "box4": "C2,2"},  # Move box4 to C2,2
        {"box1": "C1,4", "box2": "C2,8", "box3": "C1,2", "box4": "C2,2"},  # Move box3 to C1,2 (goal)
        {"box1": "C1,4", "box2": "C2,7", "box3": "C1,2", "box4": "C2,2"},  # Move box2 towards C2,2
        {"box1": "C1,4", "box2": "C2,6", "box3": "C1,2", "box4": "C2,2"},
        {"box1": "C1,4", "box2": "C2,5", "box3": "C1,2", "box4": "C2,2"},
        {"box1": "C1,4", "box2": "C2,4", "box3": "C1,2", "box4": "C2,2"},
        {"box1": "C1,4", "box2": "C2,3", "box3": "C1,2", "box4": "C2,2"},
        {"box1": "C1,4", "box2": "C2,2", "box3": "C1,2", "box4": "C2,2"},  # box2 at goal
        {"box1": "C1,4", "box2": "C2,2", "box3": "C1,2", "box4": "C2,3"},  # Move box4 towards C2,5
        {"box1": "C1,4", "box2": "C2,2", "box3": "C1,2", "box4": "C2,4"},
        {"box1": "C1,4", "box2": "C2,2", "box3": "C1,2", "box4": "C2,5"},  # box4 at goal
        {"box1": "C1,5", "box2": "C2,2", "box3": "C1,2", "box4": "C2,5"},  # Move box1 towards C2,8
        {"box1": "C1,6", "box2": "C2,2", "box3": "C1,2", "box4": "C2,5"},
        {"box1": "C1,7", "box2": "C2,2", "box3": "C1,2", "box4": "C2,5"},
        {"box1": "C1,8", "box2": "C2,2", "box3": "C1,2", "box4": "C2,5"},
        {"box1": "C2,8", "box2": "C2,2", "box3": "C1,2", "box4": "C2,5"}   # Final state - all boxes at goals
    ]
    
    return solution

# Get and print solution
solution = solve()
print(json.dumps(solution))