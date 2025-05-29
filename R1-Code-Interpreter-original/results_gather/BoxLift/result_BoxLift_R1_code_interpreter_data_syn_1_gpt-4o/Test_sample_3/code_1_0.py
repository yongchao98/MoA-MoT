from itertools import combinations

boxes = [159, 197, 211, 302, 51, 251, 299, 229, 384, 196, 377, 283, 127, 312, 228, 280, 153, 199, 169, 151]
lifters = [52, 80, 175, 75, 136, 143]

# Sort boxes in descending order
boxes.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes[:]):  # Use a copy of the list to modify it during iteration
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