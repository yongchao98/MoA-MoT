from itertools import combinations

boxes = [280, 32, 168, 120, 112, 285, 46, 33, 148, 74, 295, 279, 207, 138, 153, 63]
lifters = [53, 92, 104, 86, 106, 68]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes and len(steps) < 6:
    step = []
    used_lifters = [False] * len(lifters)
    
    for box in boxes[:]:
        # Try all combinations of lifters to lift the box
        for r in range(1, len(lifters) + 1):
            for combo in combinations(range(len(lifters)), r):
                if all(not used_lifters[i] for i in combo):
                    total_capacity = sum(lifters[i] for i in combo)
                    if total_capacity >= box:
                        step.append((box, list(combo)))
                        boxes.remove(box)
                        for i in combo:
                            used_lifters[i] = True
                        break
            if box not in boxes:
                break
    
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")