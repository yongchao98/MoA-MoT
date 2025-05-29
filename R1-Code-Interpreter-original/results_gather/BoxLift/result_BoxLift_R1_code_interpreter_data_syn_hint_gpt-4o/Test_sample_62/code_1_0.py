from itertools import combinations

boxes = [91, 207, 152, 47, 209, 54, 251, 176, 194, 221, 152, 141, 128, 159, 57, 184]
lifters = [149, 131, 113, 109, 124]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []

def can_lift(box, lifters, used):
    for r in range(1, len(lifters) + 1):
        for combo in combinations(enumerate(lifters), r):
            indices, weights = zip(*combo)
            if sum(weights) >= box and all(not used[i] for i in indices):
                return indices
    return None

while boxes and len(steps) < 6:
    step = []
    used_lifters = [False] * len(lifters)
    remaining_boxes = []
    
    for box in boxes:
        lifter_indices = can_lift(box, lifters, used_lifters)
        if lifter_indices is not None:
            step.append((box, list(lifter_indices)))
            for i in lifter_indices:
                used_lifters[i] = True
        else:
            remaining_boxes.append(box)
    
    steps.append(step)
    boxes = remaining_boxes

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")