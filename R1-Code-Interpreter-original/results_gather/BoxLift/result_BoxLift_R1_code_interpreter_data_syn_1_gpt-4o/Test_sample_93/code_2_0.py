from itertools import combinations

boxes = [64, 70, 386, 351, 113, 77, 314, 333, 266, 399, 193, 44, 181, 200, 238, 175, 370, 118, 337, 134]
lifters = [140, 115, 159, 147, 129, 112]

# Sort boxes in descending order
boxes.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    remaining_boxes = []
    
    for box in boxes:
        found = False
        for r in range(1, len(lifters) + 1):
            for combo in combinations(range(len(lifters)), r):
                if all(not used_lifters[i] for i in combo) and sum(lifters[i] for i in combo) >= box:
                    for i in combo:
                        used_lifters[i] = True
                    step.append((box, list(combo)))
                    found = True
                    break
            if found:
                break
        if not found:
            remaining_boxes.append(box)
    
    boxes = remaining_boxes
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")