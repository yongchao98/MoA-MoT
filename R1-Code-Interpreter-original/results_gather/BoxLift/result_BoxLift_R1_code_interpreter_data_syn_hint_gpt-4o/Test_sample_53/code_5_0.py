from itertools import combinations

# Box weights and lifter capacities
boxes = [108, 206, 326, 343, 45, 374, 227, 366, 110, 136, 343, 298, 201, 168, 191, 243, 365, 375, 294, 317]
lifters = [147, 80, 144, 133, 69, 200]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Function to find the minimum number of steps
def find_min_steps(boxes, lifters):
    steps = []
    while boxes:
        step = []
        remaining_lifters = lifters[:]
        used_indices = set()
        for box in boxes[:]:
            # Try to lift the box with available lifters
            best_combination = None
            for r in range(1, len(remaining_lifters) + 1):
                for combo in combinations(enumerate(remaining_lifters), r):
                    indices, selected_lifters = zip(*combo)
                    if sum(selected_lifters) >= box and not used_indices.intersection(indices):
                        if best_combination is None or len(selected_lifters) < len(best_combination[1]):
                            best_combination = (indices, selected_lifters)
            if best_combination:
                indices, selected_lifters = best_combination
                step.append((box, list(indices)))
                used_indices.update(indices)
                boxes.remove(box)
        steps.append(step)
    return steps

# Get the steps
steps = find_min_steps(boxes, lifters)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")