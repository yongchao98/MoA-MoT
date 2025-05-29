from itertools import combinations

# Box weights and lifter capacities
boxes = [41, 31, 60, 20, 11, 26, 52, 98]
lifters = [49, 55, 63]

# Function to find the best combination of lifters for a given box
def find_lifter_combination(box_weight, available_lifters):
    for r in range(1, len(available_lifters) + 1):
        for combo in combinations(enumerate(available_lifters), r):
            indices, capacities = zip(*combo)
            if sum(capacities) >= box_weight:
                return list(indices)
    return None

# Function to assign boxes to lifters in steps
def assign_boxes_to_lifters(boxes, lifters):
    steps = []
    remaining_boxes = sorted(boxes, reverse=True)
    
    while remaining_boxes:
        step = []
        used_lifters = set()
        
        for box in remaining_boxes[:]:
            available_lifters = [(i, l) for i, l in enumerate(lifters) if i not in used_lifters]
            lifter_indices = find_lifter_combination(box, [l for _, l in available_lifters])
            if lifter_indices is not None:
                step.append((box, lifter_indices))
                used_lifters.update(lifter_indices)
                remaining_boxes.remove(box)
        
        if not step:
            break
        
        steps.append(step)
    
    return steps

# Get the steps
steps = assign_boxes_to_lifters(boxes, lifters)

# Print the steps
for i, step in enumerate(steps, 1):
    print(f"Step {i}: {step}")