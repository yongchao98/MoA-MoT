from itertools import combinations

boxes = [38, 97, 21, 92, 195, 187, 95, 140, 37, 85, 25, 165]
lifters = [96, 32, 88, 48]

# Sort boxes in descending order
boxes.sort(reverse=True)

steps = []
while boxes and len(steps) < 6:
    step = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes[:]):  # Iterate over a copy of the list
        # Try to find a single lifter or a combination of lifters to lift the box
        for r in range(1, len(lifters) + 1):
            for combo in combinations(range(len(lifters)), r):
                if all(not used_lifters[j] for j in combo) and sum(lifters[j] for j in combo) >= box:
                    step.append((box, list(combo)))
                    for j in combo:
                        used_lifters[j] = True
                    boxes.remove(box)
                    break
            else:
                continue
            break
    
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")