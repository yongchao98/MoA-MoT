from itertools import combinations

# Box weights and lifter capacities
boxes = [97, 164, 187, 145, 255, 257, 240, 55, 261, 185, 252, 257, 206, 298, 121, 202]
lifters = [131, 106, 144, 121, 109, 52]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Function to find the minimum number of steps
def assign_lifters_to_boxes(boxes, lifters):
    steps = []
    while boxes:
        step = []
        used_lifters = [False] * len(lifters)
        for box in boxes[:]:
            for r in range(1, len(lifters) + 1):
                for combo in combinations(enumerate(lifters), r):
                    indices, capacities = zip(*combo)
                    if sum(capacities) >= box and all(not used_lifters[i] for i in indices):
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
steps = assign_lifters_to_boxes(boxes, lifters)

# Print the steps
for i, step in enumerate(steps, 1):
    print(f"Step {i}: {step}")