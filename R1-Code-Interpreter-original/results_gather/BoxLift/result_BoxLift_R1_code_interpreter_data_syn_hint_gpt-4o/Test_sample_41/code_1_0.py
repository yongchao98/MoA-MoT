from itertools import combinations

boxes = [88, 128, 39, 30, 80, 60, 162, 154, 191, 150, 178, 32]
lifters = [95, 111, 115, 35]

# Sort boxes in descending order
boxes.sort(reverse=True)

steps = []
remaining_boxes = boxes[:]

while remaining_boxes and len(steps) < 5:
    step = []
    used_lifters = [False] * len(lifters)
    
    # Try to lift as many boxes as possible in this step
    for box in remaining_boxes[:]:
        # Find a combination of lifters that can lift the box
        for r in range(1, len(lifters) + 1):
            for lifter_indices in combinations(range(len(lifters)), r):
                if all(not used_lifters[i] for i in lifter_indices):
                    total_capacity = sum(lifters[i] for i in lifter_indices)
                    if total_capacity >= box:
                        step.append((box, list(lifter_indices)))
                        for i in lifter_indices:
                            used_lifters[i] = True
                        remaining_boxes.remove(box)
                        break
            else:
                continue
            break
    
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")