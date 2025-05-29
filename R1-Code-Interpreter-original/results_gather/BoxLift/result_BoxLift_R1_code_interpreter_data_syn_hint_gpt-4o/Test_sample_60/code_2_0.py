from itertools import combinations

# Define the weights of the boxes and the capacities of the lifters
boxes = [21, 96, 83, 73, 100, 73, 55, 25]
lifters = [62, 79, 54]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Initialize the steps
steps = []

# Function to find the best combination of lifters for a given box
def find_lifters_for_box(box_weight, available_lifters):
    # Check all combinations of lifters
    for r in range(1, len(available_lifters) + 1):
        for combo in combinations(enumerate(available_lifters), r):
            indices, capacities = zip(*combo)
            if sum(capacities) >= box_weight:
                return list(indices)
    return []

# Lift boxes in steps
while boxes:
    step = []
    used_lifters = set()
    for box in boxes[:]:
        available_lifters = [lifters[i] for i in range(len(lifters)) if i not in used_lifters]
        lifter_indices = find_lifters_for_box(box, available_lifters)
        if lifter_indices:
            step.append((box, lifter_indices))
            used_lifters.update(lifter_indices)
            boxes.remove(box)
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")