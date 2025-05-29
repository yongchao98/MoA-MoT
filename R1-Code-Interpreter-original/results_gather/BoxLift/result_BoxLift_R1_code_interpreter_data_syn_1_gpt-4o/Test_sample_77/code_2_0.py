from itertools import combinations

boxes = [253, 162, 171, 57, 68, 142, 196, 83, 94, 287, 256, 236, 153, 79, 111, 42]
lifters = [61, 62, 122, 75, 143, 131]

# Sort boxes in descending order
boxes.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    
    for box in boxes[:]:  # Use a copy of the list to iterate
        # Try to find a combination of lifters that can lift the box
        for r in range(1, len(lifters) + 1):
            for combo in combinations(enumerate(lifters), r):
                indices, weights = zip(*combo)
                if sum(weights) >= box and all(not used_lifters[i] for i in indices):
                    step.append((box, list(indices)))
                    for i in indices:
                        used_lifters[i] = True
                    boxes.remove(box)
                    break
            else:
                continue
            break
    
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")