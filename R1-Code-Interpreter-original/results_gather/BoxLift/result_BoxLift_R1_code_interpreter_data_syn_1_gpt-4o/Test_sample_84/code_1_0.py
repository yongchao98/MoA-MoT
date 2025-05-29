from itertools import combinations

boxes = [281, 159, 239, 61, 219, 44, 225, 192, 74, 369, 318, 366, 62, 186, 373, 243, 131, 172, 66, 167, 400, 330, 351, 245]
lifters = [164, 182, 169, 127, 118, 60]

# Sort boxes in descending order
boxes.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes[:]):
        # Try to find a combination of lifters to lift the box
        lifter_indices = []
        for r in range(1, len(lifters) + 1):
            for combo in combinations(range(len(lifters)), r):
                if not any(used_lifters[j] for j in combo):
                    if sum(lifters[j] for j in combo) >= box:
                        lifter_indices = list(combo)
                        break
            if lifter_indices:
                break
        
        if lifter_indices:
            step.append((box, lifter_indices))
            boxes.remove(box)
            for j in lifter_indices:
                used_lifters[j] = True
    
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")