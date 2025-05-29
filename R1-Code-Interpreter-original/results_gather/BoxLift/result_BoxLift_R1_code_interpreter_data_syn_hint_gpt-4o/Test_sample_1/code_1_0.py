from itertools import combinations

boxes = [71, 58, 105, 85, 123, 125, 84, 105, 141, 73, 102, 46]
lifters = [77, 45, 120, 109]

# Sort boxes in descending order
boxes.sort(reverse=True)

steps = []
while boxes and len(steps) < 5:
    step = []
    remaining_lifters = lifters[:]
    used_lifters = set()
    
    for box in boxes[:]:
        # Try to find a single lifter or a combination of lifters to lift the box
        found = False
        for r in range(1, len(remaining_lifters) + 1):
            for combo in combinations(enumerate(remaining_lifters), r):
                indices, lifter_combo = zip(*combo)
                if sum(lifter_combo) >= box and not any(i in used_lifters for i in indices):
                    step.append((box, list(indices)))
                    for i in indices:
                        used_lifters.add(i)
                    boxes.remove(box)
                    found = True
                    break
            if found:
                break
    
    steps.append(step)

# Print the steps
for idx, step in enumerate(steps):
    print(f"Step {idx + 1}: {step}")