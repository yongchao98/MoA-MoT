from itertools import combinations

# Box weights and lifter capacities
boxes = [181, 52, 82, 130, 92, 195, 101, 106, 189, 42, 196, 46]
lifters = [45, 30, 52, 104, 30]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Function to find the minimum number of steps
def find_min_steps(boxes, lifters):
    steps = []
    while boxes:
        step = []
        used_lifters = [False] * len(lifters)
        for box in boxes[:]:
            # Try to lift the box with available lifters
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
    return steps

# Get the steps
steps = find_min_steps(boxes, lifters)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")