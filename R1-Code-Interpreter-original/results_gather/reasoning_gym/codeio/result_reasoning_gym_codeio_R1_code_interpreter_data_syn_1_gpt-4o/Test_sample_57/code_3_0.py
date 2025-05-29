import numpy as np

def find_feasible_input():
    # Given output
    output = [8.188638861262342, -3.1143668514063894, 11.308441967021771, 12.098966938050486, 9.984426035934856]
    
    # Assume B and a
    B = 0.5
    a = 0.1
    
    # Calculate constraints and global best
    constraints = []
    global_best = []
    
    for i, out in enumerate(output):
        # Assume L = 10 for simplicity
        L = 10
        # Calculate global best assuming e = 0
        gb = out / (B * L)
        global_best.append(gb)
        # Calculate constraints assuming a range around the global best
        constraints.append([gb - 5, gb + 5])
    
    return {
        "constraints": constraints,
        "global_best": global_best,
        "B": B,
        "a": a
    }

print(find_feasible_input())