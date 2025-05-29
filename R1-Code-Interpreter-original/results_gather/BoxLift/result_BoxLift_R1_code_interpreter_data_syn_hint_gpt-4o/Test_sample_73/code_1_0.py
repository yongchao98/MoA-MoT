from itertools import combinations

boxes = [75, 210, 297, 242, 276, 134, 383, 358, 146, 394, 59, 350, 94, 152, 250, 136, 116, 342, 283, 87]
lifters = [101, 62, 125, 126, 138, 54]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    
    for box in boxes[:]:  # Iterate over a copy of the list
        # Try to find a combination of lifters to lift the box
        lifter_indices = []
        for r in range(1, len(lifters) + 1):
            for combo in combinations(range(len(lifters)), r):
                if all(not used_lifters[i] for i in combo):
                    if sum(lifters[i] for i in combo) >= box:
                        lifter_indices = list(combo)
                        break
            if lifter_indices:
                break
        
        if lifter_indices:
            step.append((box, lifter_indices))
            for i in lifter_indices:
                used_lifters[i] = True
            boxes.remove(box)
    
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")