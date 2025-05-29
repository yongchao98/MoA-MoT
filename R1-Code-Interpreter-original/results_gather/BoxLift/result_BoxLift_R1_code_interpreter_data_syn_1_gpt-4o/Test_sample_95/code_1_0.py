from itertools import combinations

# Box weights and lifter capacities
boxes = [79, 48, 16, 95, 67, 41, 62, 22]
lifters = [52, 41, 41, 67]

# Sort boxes in descending order to prioritize heavier boxes
boxes.sort(reverse=True)

# Function to find the best combination of lifters for a given box
def find_lifters_for_box(box_weight, available_lifters):
    for r in range(1, len(available_lifters) + 1):
        for combo in combinations(available_lifters, r):
            if sum(combo) >= box_weight:
                return combo
    return None

# Function to assign lifters to boxes in steps
def assign_lifters_to_boxes(boxes, lifters):
    steps = []
    while boxes:
        step = []
        available_lifters = lifters[:]
        used_indices = set()
        for box in boxes[:]:
            lifter_combo = find_lifters_for_box(box, available_lifters)
            if lifter_combo:
                lifter_indices = []
                for l in lifter_combo:
                    index = available_lifters.index(l)
                    while index in used_indices:
                        index = available_lifters.index(l, index + 1)
                    lifter_indices.append(index)
                    used_indices.add(index)
                step.append((box, lifter_indices))
                boxes.remove(box)
                for index in lifter_indices:
                    available_lifters[index] = 0  # Mark as used
        steps.append(step)
        if len(steps) > 3:
            return "Not possible to lift all boxes in 3 steps"
    return steps

# Get the steps
steps = assign_lifters_to_boxes(boxes, lifters)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")