from itertools import combinations

boxes = [351, 64, 363, 78, 283, 68, 211, 200, 288, 295, 319, 271, 305, 219, 310, 365, 390, 56, 396, 332, 130, 136, 221, 386]
lifters = [134, 184, 163, 144, 109, 179, 64]

# Sort boxes in descending order
boxes.sort(reverse=True)

steps = []
while boxes and len(steps) < 8:
    step = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes[:]):
        # Try to find a combination of lifters that can lift the box
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

for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")