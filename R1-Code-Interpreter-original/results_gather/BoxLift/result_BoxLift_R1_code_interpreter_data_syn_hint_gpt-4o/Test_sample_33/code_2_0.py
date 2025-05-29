from itertools import combinations

# Box weights and lifter capacities
boxes = [88, 47, 209, 393, 59, 342, 259, 103, 304, 239, 141, 90, 256, 176, 297, 274, 154, 159, 262, 183]
lifters = [173, 188, 50, 119, 177, 86]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Function to find the minimum number of steps
def assign_lifters_to_boxes(boxes, lifters):
    steps = []
    while boxes:
        step = []
        used_lifters = [False] * len(lifters)
        for i, box in enumerate(boxes[:]):
            # Try to lift the box with available lifters
            lifter_indices = []
            remaining_weight = box
            for r in range(1, len(lifters) + 1):
                for combo in combinations(range(len(lifters)), r):
                    if all(not used_lifters[j] for j in combo) and sum(lifters[j] for j in combo) >= remaining_weight:
                        lifter_indices = list(combo)
                        for j in combo:
                            used_lifters[j] = True
                        break
                if lifter_indices:
                    break
            if lifter_indices:
                step.append((box, lifter_indices))
                boxes.remove(box)
        steps.append(step)
        if len(steps) > 8:
            return "Not possible to lift all boxes in 8 steps"
    return steps

# Get the steps
steps = assign_lifters_to_boxes(boxes, lifters)

# Print the steps
if isinstance(steps, str):
    print(steps)
else:
    for i, step in enumerate(steps):
        print(f"Step {i+1}: {step}")