from itertools import combinations

boxes = [99, 159, 55, 86, 215, 56, 271, 167, 212, 87, 78, 211, 246, 153, 227, 287]
lifters = [143, 50, 70, 79, 136]

# Sort boxes in descending order
boxes.sort(reverse=True)

steps = []
while boxes and len(steps) < 7:
    step = []
    remaining_lifters = lifters[:]
    used_lifters = set()
    
    for box in boxes[:]:
        # Try to find a combination of lifters to lift the box
        for r in range(1, len(remaining_lifters) + 1):
            for combo in combinations(enumerate(remaining_lifters), r):
                indices, capacities = zip(*combo)
                if sum(capacities) >= box and not used_lifters.intersection(indices):
                    step.append((box, list(indices)))
                    boxes.remove(box)
                    used_lifters.update(indices)
                    break
            if box not in boxes:
                break
    
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")