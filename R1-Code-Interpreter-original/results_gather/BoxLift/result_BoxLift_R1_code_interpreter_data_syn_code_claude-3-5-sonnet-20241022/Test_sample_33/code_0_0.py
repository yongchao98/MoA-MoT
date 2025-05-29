import itertools

# Input data
boxes = [88, 47, 209, 393, 59, 342, 259, 103, 304, 239, 141, 90, 256, 176, 297, 274, 154, 159, 262, 183]
lifters = [173, 188, 50, 119, 177, 86]

def find_lifter_combinations(box_weight, available_lifters, lifter_indices):
    # Try all possible combinations of lifters
    for r in range(1, len(available_lifters) + 1):
        for comb in itertools.combinations(range(len(available_lifters)), r):
            total_capacity = sum(available_lifters[i] for i in comb)
            if total_capacity >= box_weight:
                # Convert to original lifter indices
                return [lifter_indices[i] for i in comb]
    return None

def solve_box_lifting():
    # Sort boxes in descending order
    sorted_boxes = sorted(enumerate(boxes), key=lambda x: x[1], reverse=True)
    steps = []
    remaining_boxes = list(sorted_boxes)
    
    while remaining_boxes:
        step = []
        available_lifters = list(enumerate(lifters))
        used_lifters = set()
        
        # Try to assign boxes to lifters
        i = 0
        while i < len(remaining_boxes) and available_lifters:
            box_idx, box_weight = remaining_boxes[i]
            
            # Get available lifter indices
            available_indices = [idx for idx, _ in available_lifters]
            available_capacities = [cap for _, cap in available_lifters]
            
            # Find lifters for this box
            lifter_combination = find_lifter_combinations(box_weight, available_capacities, available_indices)
            
            if lifter_combination:
                # Add assignment to step
                step.append((box_weight, lifter_combination))
                # Remove used lifters
                available_lifters = [(idx, cap) for idx, cap in available_lifters 
                                   if idx not in lifter_combination]
                # Remove box from remaining boxes
                remaining_boxes.pop(i)
            else:
                i += 1
        
        steps.append(step)
        
        if len(steps) > 8:
            return None
    
    return steps

# Solve and format output
solution = solve_box_lifting()
if solution:
    output = ""
    for i, step in enumerate(solution, 1):
        output += f"Step {i}: {step}\n"
    print(f"<<{output}>>")
else:
    print("No solution found within 8 steps.")