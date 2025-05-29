from itertools import combinations

# Box weights and lifter capacities
boxes = [58, 45, 94, 26, 49, 153, 61, 23, 170, 143, 75, 166]
lifters = [57, 61, 104, 98, 70]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Function to find the optimal assignment of lifters to boxes
def assign_lifters_to_boxes(boxes, lifters):
    steps = []
    while boxes:
        step = []
        used_lifters = [False] * len(lifters)
        remaining_lifters = [(i, lifter) for i, lifter in enumerate(lifters)]
        
        for box in boxes[:]:
            # Try to find a combination of lifters to lift the box
            for r in range(1, len(lifters) + 1):
                for combo in combinations(remaining_lifters, r):
                    indices, capacities = zip(*combo)
                    if sum(capacities) >= box and all(not used_lifters[i] for i in indices):
                        step.append((box, list(indices)))
                        for i in indices:
                            used_lifters[i] = True
                        boxes.remove(box)
                        break
                else:
                    continue
                break
        steps.append(step)
        if len(steps) > 4:
            break
    return steps

# Get the steps
steps = assign_lifters_to_boxes(boxes, lifters)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")