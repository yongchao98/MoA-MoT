from itertools import combinations

boxes = [356, 292, 240, 155, 142, 223, 178, 231, 254, 315, 236, 328, 203, 397, 44, 200, 80, 131, 175, 216]
lifters = [155, 116, 53, 98, 96, 156, 167]

# Sort boxes in descending order
boxes.sort(reverse=True)

steps = []
while boxes and len(steps) < 7:
    step = []
    used_lifters = [False] * len(lifters)
    
    for box in boxes[:]:  # Use a copy of the list to modify it during iteration
        # Try to find a combination of lifters that can lift the box
        for r in range(1, len(lifters) + 1):
            for combo in combinations(range(len(lifters)), r):
                if all(not used_lifters[i] for i in combo):
                    total_capacity = sum(lifters[i] for i in combo)
                    if total_capacity >= box:
                        step.append((box, list(combo)))
                        for i in combo:
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