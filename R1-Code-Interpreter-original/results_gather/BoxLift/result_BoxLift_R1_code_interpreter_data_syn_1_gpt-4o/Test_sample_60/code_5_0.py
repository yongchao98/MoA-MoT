from itertools import combinations

# Box weights and lifter capacities
boxes = [21, 96, 83, 73, 100, 73, 55, 25]
lifters = [62, 79, 54]

# Function to find the best combination of lifters to lift a box
def find_lifter_combination(box_weight, available_lifters):
    for r in range(1, len(available_lifters) + 1):
        for combo in combinations(enumerate(available_lifters), r):
            indices, capacities = zip(*combo)
            if sum(capacities) >= box_weight:
                return list(indices)
    return None

# Function to assign lifters to boxes in steps
def assign_lifters_to_boxes(boxes, lifters):
    steps = []
    remaining_boxes = sorted(boxes, reverse=True)
    
    while remaining_boxes:
        step = []
        used_lifters = set()
        
        for box_weight in remaining_boxes[:]:
            available_lifters = [l for i, l in enumerate(lifters) if i not in used_lifters]
            lifter_indices = find_lifter_combination(box_weight, available_lifters)
            if lifter_indices is not None:
                step.append((box_weight, lifter_indices))
                used_lifters.update(lifter_indices)
                remaining_boxes.remove(box_weight)
        
        steps.append(step)
    
    return steps

# Get the steps
steps = assign_lifters_to_boxes(boxes, lifters)

# Print the steps
for i, step in enumerate(steps, 1):
    print(f"Step {i}: {step}")