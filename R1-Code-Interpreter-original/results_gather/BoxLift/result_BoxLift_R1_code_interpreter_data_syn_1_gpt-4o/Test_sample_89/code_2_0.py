from itertools import combinations

boxes = [227, 106, 388, 173, 63, 178, 261, 205, 323, 124, 386, 81, 71, 127, 134, 212, 150, 114, 41, 277, 123, 152, 47, 47]
lifters = [166, 106, 122, 181, 102, 119, 147]

# Sort boxes in descending order
boxes.sort(reverse=True)

steps = []
while boxes and len(steps) < 6:
    step = []
    used_lifters = [False] * len(lifters)
    remaining_boxes = []
    
    for box in boxes:
        # Try to find a combination of lifters to lift the box
        found = False
        for r in range(1, len(lifters) + 1):
            for combo in combinations(range(len(lifters)), r):
                if not any(used_lifters[i] for i in combo) and sum(lifters[i] for i in combo) >= box:
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