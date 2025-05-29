from itertools import combinations

boxes = [21, 25, 138, 168, 105, 127, 154, 55, 155, 125, 121, 58]
lifters = [78, 98, 51, 75, 43]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
    step = []
    remaining_lifters = lifters[:]
    used_lifters = [False] * len(lifters)
    for box in boxes[:]:
        # Try to find a combination of lifters to lift the box
        found = False
        for r in range(1, len(remaining_lifters) + 1):
            for combo in combinations(enumerate(remaining_lifters), r):
                indices, selected_lifters = zip(*combo)
                if sum(selected_lifters) >= box and all(not used_lifters[i] for i in indices):
                    step.append((box, list(indices)))
                    for i in indices:
                        used_lifters[i] = True
                    boxes.remove(box)
                    found = True
                    break
            if found:
                break
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")