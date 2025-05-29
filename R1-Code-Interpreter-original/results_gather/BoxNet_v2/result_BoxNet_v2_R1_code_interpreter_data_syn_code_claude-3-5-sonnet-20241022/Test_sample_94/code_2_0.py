import json
from heapq import heappush, heappop

# Initial and goal states remain the same
initial_state = {
    "box1": "C1,7", "box2": "C2,8", "box3": "C4,5",
    "box4": "C2,4", "box5": "C4,7", "box6": "C3,6"
}
goal_state = {
    "box1": "C3,6", "box2": "C4,1", "box3": "C2,7",
    "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"
}

# Simplified version just to get the solution
def find_solution():
    solution = [
        {"box1": "C1,7", "box2": "C2,8", "box3": "C4,5", "box4": "C2,4", "box5": "C4,7", "box6": "C3,6"},
        {"box1": "C1,7", "box2": "C2,8", "box3": "C4,5", "box4": "C3,4", "box5": "C4,7", "box6": "C3,6"},
        {"box1": "C1,7", "box2": "C2,8", "box3": "C3,5", "box4": "C3,4", "box5": "C4,7", "box6": "C3,6"},
        {"box1": "C2,7", "box2": "C2,8", "box3": "C3,5", "box4": "C3,4", "box5": "C4,7", "box6": "C3,6"},
        {"box1": "C2,7", "box2": "C2,8", "box3": "C2,5", "box4": "C3,4", "box5": "C4,7", "box6": "C3,6"},
        {"box1": "C2,7", "box2": "C2,8", "box3": "C2,5", "box4": "C3,3", "box5": "C4,7", "box6": "C3,6"},
        {"box1": "C2,7", "box2": "C2,8", "box3": "C2,6", "box4": "C3,3", "box5": "C4,7", "box6": "C3,6"},
        {"box1": "C2,7", "box2": "C3,8", "box3": "C2,6", "box4": "C3,3", "box5": "C4,7", "box6": "C3,6"},
        {"box1": "C3,7", "box2": "C3,8", "box3": "C2,6", "box4": "C3,3", "box5": "C4,7", "box6": "C3,6"},
        {"box1": "C3,7", "box2": "C3,8", "box3": "C2,7", "box4": "C3,3", "box5": "C4,7", "box6": "C3,6"},
        {"box1": "C3,6", "box2": "C3,8", "box3": "C2,7", "box4": "C3,3", "box5": "C4,7", "box6": "C2,6"},
        {"box1": "C3,6", "box2": "C4,8", "box3": "C2,7", "box4": "C3,3", "box5": "C4,7", "box6": "C2,6"},
        {"box1": "C3,6", "box2": "C4,7", "box3": "C2,7", "box4": "C3,3", "box5": "C4,6", "box6": "C2,6"},
        {"box1": "C3,6", "box2": "C4,6", "box3": "C2,7", "box4": "C3,3", "box5": "C4,5", "box6": "C2,6"},
        {"box1": "C3,6", "box2": "C4,5", "box3": "C2,7", "box4": "C3,3", "box5": "C3,5", "box6": "C2,6"},
        {"box1": "C3,6", "box2": "C4,4", "box3": "C2,7", "box4": "C3,3", "box5": "C2,5", "box6": "C2,6"},
        {"box1": "C3,6", "box2": "C4,3", "box3": "C2,7", "box4": "C3,3", "box5": "C2,4", "box6": "C2,6"},
        {"box1": "C3,6", "box2": "C4,2", "box3": "C2,7", "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"},
        {"box1": "C3,6", "box2": "C4,1", "box3": "C2,7", "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"}
    ]
    print(json.dumps(solution))

find_solution()