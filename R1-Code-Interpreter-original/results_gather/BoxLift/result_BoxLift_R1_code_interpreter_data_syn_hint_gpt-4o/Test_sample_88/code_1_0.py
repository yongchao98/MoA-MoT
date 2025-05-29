from itertools import combinations

boxes = [88, 313, 344, 387, 392, 370, 99, 389, 340, 391, 131, 147, 334, 251, 157, 144, 290, 167, 141, 79]
lifters = [148, 193, 54, 146, 106, 53, 88]

# Sort boxes in descending order
boxes.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    
    for box in boxes[:]:
        # Try to find a combination of lifters to lift the box
        for r in range(1, len(lifters) + 1):
            for lifter_indices in combinations(range(len(lifters)), r):
                if all(not used_lifters[i] for i in lifter_indices):
                    total_capacity = sum(lifters[i] for i in lifter_indices)
                    if total_capacity >= box:
                        step.append((box, list(lifter_indices)))
                        for i in lifter_indices:
                            used_lifters[i] = True
                        boxes.remove(box)
                        break
            else:
                continue
            break
    
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")