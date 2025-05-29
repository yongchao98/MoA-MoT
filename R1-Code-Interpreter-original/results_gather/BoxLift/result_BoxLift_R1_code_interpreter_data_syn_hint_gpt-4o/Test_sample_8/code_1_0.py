from itertools import combinations

# Box weights and lifter capacities
boxes = [362, 202, 192, 149, 357, 335, 209, 193, 334, 375, 374, 167, 276, 70, 282, 201, 242, 158, 236, 212]
lifters = [174, 130, 182, 186, 80, 149, 87]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

# Function to find the optimal assignment of lifters to boxes
def assign_lifters_to_boxes(boxes, lifters):
    steps = []
    while boxes:
        step = []
        used_lifters = [False] * len(lifters)
        for i, box in enumerate(boxes[:]):
            # Try to find a combination of lifters to lift the box
            lifter_indices = []
            for r in range(1, len(lifters) + 1):
                for combo in combinations(range(len(lifters)), r):
                    if not any(used_lifters[j] for j in combo):
                        if sum(lifters[j] for j in combo) >= box:
                            lifter_indices = list(combo)
                            break
                if lifter_indices:
                    break
            if lifter_indices:
                step.append((box, lifter_indices))
                boxes.remove(box)
                for j in lifter_indices:
                    used_lifters[j] = True
        steps.append(step)
    return steps

# Get the steps
steps = assign_lifters_to_boxes(boxes, lifters)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")