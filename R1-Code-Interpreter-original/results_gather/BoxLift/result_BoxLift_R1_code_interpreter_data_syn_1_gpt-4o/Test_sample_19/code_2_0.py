from itertools import combinations

boxes = [55, 77, 274, 297, 70, 312, 108, 84, 340, 273, 157, 386, 95, 223, 311, 306, 287, 91, 220, 296, 320, 237, 97, 73]
lifters = [167, 138, 133, 87, 142, 167, 140]

# Sort boxes in descending order
boxes.sort(reverse=True)

steps = []
while boxes and len(steps) < 7:
    step = []
    used_lifters = [False] * len(lifters)
    remaining_boxes = []
    
    for box in boxes:
        # Try to find a combination of lifters to lift the box
        found = False
        for r in range(1, len(lifters) + 1):
            for combo in combinations(range(len(lifters)), r):
                if not any(used_lifters[i] for i in combo):
                    if sum(lifters[i] for i in combo) >= box:
                        step.append((box, list(combo)))
                        for i in combo:
                            used_lifters[i] = True
                        found = True
                        break
            if found:
                break
        if not found:
            remaining_boxes.append(box)
    
    steps.append(step)
    boxes = remaining_boxes

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")