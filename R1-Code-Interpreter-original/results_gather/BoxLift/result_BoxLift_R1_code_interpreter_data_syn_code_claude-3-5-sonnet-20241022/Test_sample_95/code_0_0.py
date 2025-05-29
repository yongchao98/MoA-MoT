from itertools import combinations

def can_lift(lifters, weight):
    # Try all possible combinations of lifters to see if they can lift the weight
    for i in range(1, len(lifters) + 1):
        for combo in combinations(range(len(lifters)), i):
            if sum(lifters[j] for j in combo) >= weight:
                return list(combo)
    return None

def find_solution():
    boxes = [79, 48, 16, 95, 67, 41, 62, 22]
    lifter_capacities = [52, 41, 41, 67]
    boxes_with_indices = list(enumerate(boxes))
    solution = []
    
    while boxes_with_indices:
        step = []
        available_lifters = set(range(len(lifter_capacities)))
        
        # Sort boxes by weight (descending) to handle heaviest boxes first
        boxes_with_indices.sort(key=lambda x: x[1], reverse=True)
        
        i = 0
        while i < len(boxes_with_indices) and available_lifters:
            box_idx, weight = boxes_with_indices[i]
            
            # Try to find lifters for this box
            possible_lifters = [idx for idx in available_lifters]
            lifter_subset = can_lift([lifter_capacities[j] for j in possible_lifters], weight)
            
            if lifter_subset is not None:
                # Convert subset indices to actual lifter indices
                actual_lifters = [possible_lifters[j] for j in lifter_subset]
                step.append((weight, actual_lifters))
                available_lifters -= set(actual_lifters)
                boxes_with_indices.pop(i)
            else:
                i += 1
        
        solution.append(step)
        if len(solution) > 3:
            return None
    
    return solution

# Find and format solution
solution = find_solution()
if solution:
    formatted_solution = ""
    for i, step in enumerate(solution, 1):
        formatted_solution += f"Step {i}: {step}\n"
    print(formatted_solution.strip())
else:
    print("No solution found within 3 steps")