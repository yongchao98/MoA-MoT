from itertools import combinations

boxes = [71, 55, 11, 10, 71, 91, 14, 69]
lifters = [46, 54, 47]

# Sort boxes in descending order
boxes.sort(reverse=True)

steps = []
while boxes and len(steps) < 5:
    step = []
    used_lifters = [False] * len(lifters)
    remaining_boxes = []
    
    for box in boxes:
        # Try to find a combination of lifters to lift the box
        found = False
        for r in range(1, len(lifters) + 1):
            for lifter_indices in combinations(range(len(lifters)), r):
                if all(not used_lifters[i] for i in lifter_indices):
                    if sum(lifters[i] for i in lifter_indices) >= box:
                        step.append((box, list(lifter_indices)))
                        for i in lifter_indices:
                            used_lifters[i] = True
                        found = True
                        break
            if found:
                break
        if not found:
            remaining_boxes.append(box)
    
    steps.append(step)
    boxes = remaining_boxes

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")