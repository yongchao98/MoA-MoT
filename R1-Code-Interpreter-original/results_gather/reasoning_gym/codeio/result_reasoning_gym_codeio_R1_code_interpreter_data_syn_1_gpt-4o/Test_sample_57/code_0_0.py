import numpy as np

def find_feasible_input():
    # Given output
    output_positions = [8.188638861262342, -3.1143668514063894, 11.308441967021771, 12.098966938050486, 9.984426035934856]
    
    # Assume constraints for each dimension
    constraints = [[-10, 10] for _ in range(len(output_positions))]
    
    # Assume global best positions close to the output
    global_best = [8.0, -3.0, 11.0, 12.0, 10.0]
    
    # Typical values for B and a
    B = 0.5
    a = 0.1
    
    # Return the assumed inputs
    return {
        "constraints": constraints,
        "global_best": global_best,
        "B": B,
        "a": a
    }

# Execute the function to find feasible input
feasible_input = find_feasible_input()
print(feasible_input)