from itertools import combinations

boxes = [44, 70, 139, 184, 194, 67, 188, 33, 149, 194, 132, 117]
lifters = [110, 109, 30, 34]

# Sort boxes in descending order
boxes.sort(reverse=True)

steps = []
while boxes:
    step = []
    remaining_lifters = lifters[:]
    used_lifters = [False] * len(lifters)
    
    for box in boxes[:]:
        # Try to find a combination of lifters to lift the box
        for r in range(1, len(lifters) + 1):
            for combo in combinations(enumerate(remaining_lifters), r):
                indices, lifter_combo = zip(*combo)
                if sum(lifter_combo) >= box and all(not used_lifters[i] for i in indices):
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