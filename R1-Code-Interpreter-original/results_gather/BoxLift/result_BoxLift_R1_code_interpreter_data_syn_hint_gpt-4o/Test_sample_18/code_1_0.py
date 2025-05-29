from itertools import combinations

boxes = [161, 400, 306, 233, 179, 305, 284, 214, 288, 132, 387, 198, 274, 40, 392, 183, 275, 259, 67, 177]
lifters = [122, 116, 51, 165, 95, 150, 179]

# Sort boxes in descending order
boxes.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes):
        if not boxes:
            break
        # Try to find a combination of lifters that can lift the box
        for r in range(1, len(lifters) + 1):
            for lifter_indices in combinations(range(len(lifters)), r):
                if all(not used_lifters[idx] for idx in lifter_indices):
                    total_capacity = sum(lifters[idx] for idx in lifter_indices)
                    if total_capacity >= box:
                        for idx in lifter_indices:
                            used_lifters[idx] = True
                        step.append((box, list(lifter_indices)))
                        boxes[i] = None
                        break
            if boxes[i] is None:
                break
    
    # Remove lifted boxes
    boxes = [box for box in boxes if box is not None]
    steps.append(step)

# Print the steps
for idx, step in enumerate(steps):
    print(f"Step {idx + 1}: {step}")