from itertools import combinations

boxes = [187, 170, 257, 238, 229, 202, 239, 125, 197, 243, 110, 50, 48, 168, 42, 292]
lifters = [115, 55, 142, 79, 63]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes and len(steps) < 10:
    step = []
    used_lifters = [False] * len(lifters)
    
    for box in boxes[:]:  # Iterate over a copy of the list
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