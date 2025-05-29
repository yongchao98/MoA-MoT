from itertools import combinations

boxes = [85, 162, 147, 83, 142, 96, 200, 172, 151, 77, 59, 39]
lifters = [53, 54, 76, 48, 97]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
    step = []
    remaining_lifters = lifters[:]
    used_indices = set()
    for box in boxes[:]:
        # Try to find a combination of lifters to lift the box
        found = False
        for r in range(1, len(remaining_lifters) + 1):
            for combo in combinations(range(len(remaining_lifters)), r):
                if sum(remaining_lifters[i] for i in combo) >= box and not any(i in used_indices for i in combo):
                    # Box can be lifted, add to step
                    step.append((box, list(combo)))
                    # Mark lifters as used
                    used_indices.update(combo)
                    # Remove the box from the list
                    boxes.remove(box)
                    found = True
                    break
            if found:
                break
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")