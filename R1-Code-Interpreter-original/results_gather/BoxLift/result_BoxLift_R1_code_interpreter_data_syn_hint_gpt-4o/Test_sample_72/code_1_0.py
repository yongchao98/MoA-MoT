from itertools import combinations

boxes = [280, 109, 195, 129, 54, 216, 157, 193, 44, 233, 171, 218, 204, 120, 131, 151]
lifters = [49, 131, 66, 124, 91]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes[:]):  # Iterate over a copy of the list
        # Try to find a combination of lifters to lift the box
        for r in range(1, len(lifters) + 1):
            for lifter_indices in combinations(range(len(lifters)), r):
                if all(not used_lifters[j] for j in lifter_indices):
                    total_capacity = sum(lifters[j] for j in lifter_indices)
                    if total_capacity >= box:
                        step.append((box, list(lifter_indices)))
                        for j in lifter_indices:
                            used_lifters[j] = True
                        boxes.remove(box)
                        break
            else:
                continue
            break
    
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")