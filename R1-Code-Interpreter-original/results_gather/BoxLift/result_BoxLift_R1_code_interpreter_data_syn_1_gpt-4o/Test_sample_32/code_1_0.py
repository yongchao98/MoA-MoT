from itertools import combinations

boxes = [257, 47, 186, 194, 237, 135, 72, 204, 268, 83, 242, 41, 214, 157, 115, 242]
lifters = [77, 41, 40, 54, 132]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes):
        if not boxes:
            break
        # Try to find a combination of lifters to lift the box
        lifter_indices = []
        remaining_weight = box
        
        # Check all combinations of lifters
        for r in range(1, len(lifters) + 1):
            for combo in combinations(range(len(lifters)), r):
                if all(not used_lifters[j] for j in combo):
                    total_capacity = sum(lifters[j] for j in combo)
                    if total_capacity >= remaining_weight:
                        lifter_indices = list(combo)
                        for j in combo:
                            used_lifters[j] = True
                        break
            if lifter_indices:
                break
        
        if lifter_indices:
            step.append((box, lifter_indices))
            boxes[i] = None  # Mark this box as lifted
    
    # Remove lifted boxes
    boxes = [box for box in boxes if box is not None]
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")