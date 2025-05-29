from itertools import combinations

boxes = [161, 400, 306, 233, 179, 305, 284, 214, 288, 132, 387, 198, 274, 40, 392, 183, 275, 259, 67, 177]
lifters = [122, 116, 51, 165, 95, 150, 179]

# Sort boxes in descending order
boxes.sort(reverse=True)

steps = []
while boxes and len(steps) < 7:
    step = []
    used_lifters = [False] * len(lifters)
    remaining_boxes = []
    
    for box in boxes:
        # Try to find a combination of lifters to lift the current box
        found = False
        for r in range(1, len(lifters) + 1):
            for combo in combinations(enumerate(lifters), r):
                indices, weights = zip(*combo)
                if not any(used_lifters[i] for i in indices) and sum(weights) >= box:
                    step.append((box, list(indices)))
                    for i in indices:
                        used_lifters[i] = True
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