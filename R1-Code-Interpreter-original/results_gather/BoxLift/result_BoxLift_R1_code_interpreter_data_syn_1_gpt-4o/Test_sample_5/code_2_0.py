from itertools import combinations

# Box weights and lifter capacities
boxes = [38, 74, 55, 11, 82, 43, 34, 20]
lifters = [76, 78, 46, 63]

# Sort boxes in descending order to prioritize heavier boxes
boxes.sort(reverse=True)

# Sort lifters in descending order to prioritize stronger lifters
lifters.sort(reverse=True)

# Function to find the best combination of lifters for a given box
def find_lifters_for_box(box_weight, available_lifters):
    # Try to find a single lifter who can lift the box
    for i, lifter in enumerate(available_lifters):
        if lifter >= box_weight:
            return [i]
    
    # If no single lifter can lift the box, try combinations
    for r in range(2, len(available_lifters) + 1):
        for indices in combinations(range(len(available_lifters)), r):
            if sum(available_lifters[i] for i in indices) >= box_weight:
                return list(indices)
    
    # If no combination can lift the box, return an empty list
    return []

# Function to assign lifters to boxes in steps
def assign_lifters_to_boxes(boxes, lifters):
    steps = []
    while boxes:
        step = []
        available_lifters = lifters[:]
        used_lifters = set()
        for box in boxes[:]:
            lifter_indices = find_lifters_for_box(box, available_lifters)
            if lifter_indices and not any(i in used_lifters for i in lifter_indices):
                step.append((box, lifter_indices))
                boxes.remove(box)
                used_lifters.update(lifter_indices)
                for index in sorted(lifter_indices, reverse=True):
                    available_lifters.pop(index)
        steps.append(step)
        if len(steps) > 3:
            break
    return steps

# Get the steps for lifting all boxes
steps = assign_lifters_to_boxes(boxes, lifters)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")