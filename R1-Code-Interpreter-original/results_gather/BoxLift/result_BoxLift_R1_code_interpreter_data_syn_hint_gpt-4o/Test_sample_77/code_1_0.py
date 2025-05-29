from itertools import combinations

boxes = [253, 162, 171, 57, 68, 142, 196, 83, 94, 287, 256, 236, 153, 79, 111, 42]
lifters = [61, 62, 122, 75, 143, 131]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes and len(steps) < 6:
    step = []
    remaining_lifters = lifters[:]
    used_lifters = set()
    for box in boxes[:]:
        # Try to find a combination of lifters to lift the box
        for r in range(1, len(remaining_lifters) + 1):
            for combo in combinations(range(len(remaining_lifters)), r):
                if combo not in used_lifters:
                    total_capacity = sum(remaining_lifters[i] for i in combo)
                    if total_capacity >= box:
                        # Assign this combination to lift the box
                        step.append((box, list(combo)))
                        used_lifters.update(combo)
                        # Remove the box from the list
                        boxes.remove(box)
                        break
            else:
                continue
            break
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")