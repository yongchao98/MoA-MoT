from itertools import combinations

boxes = [151, 393, 143, 179, 88, 173, 344, 212, 308, 258, 113, 261, 371, 335, 236, 258, 124, 262, 108, 209, 232, 128, 266, 232]
lifters = [84, 94, 87, 108, 198, 95, 121]

# Sort boxes in descending order
boxes.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    remaining_boxes = []
    
    for box in boxes:
        # Try to find the smallest combination of lifters that can lift the box
        found = False
        for r in range(1, len(lifters) + 1):
            for combo in combinations(range(len(lifters)), r):
                if not any(used_lifters[i] for i in combo):
                    if sum(lifters[i] for i in combo) >= box:
                        for i in combo:
                            used_lifters[i] = True
                        step.append((box, list(combo)))
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
    print(f"Step {i+1}: {step}")