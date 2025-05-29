from itertools import combinations

boxes = [116, 246, 296, 369, 78, 275, 77, 383, 71, 155, 352, 94, 63, 168, 350, 79, 59, 252, 88, 278, 188, 383, 240, 308]
lifters = [134, 132, 150, 137, 62, 132]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes):
        if box is None:
            continue
        # Try to find a combination of lifters that can lift the box
        for r in range(1, len(lifters) + 1):
            for combo in combinations(range(len(lifters)), r):
                if all(not used_lifters[j] for j in combo):
                    total_capacity = sum(lifters[j] for j in combo)
                    if total_capacity >= box:
                        step.append((box, list(combo)))
                        for j in combo:
                            used_lifters[j] = True
                        boxes[i] = None
                        break
            if boxes[i] is None:
                break
    
    boxes = [box for box in boxes if box is not None]
    steps.append(step)

for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")