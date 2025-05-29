from itertools import combinations

boxes = [94, 177, 88, 131, 87, 113, 148, 86, 141, 83, 86, 73]
lifters = [55, 92, 35, 48, 96]

# Sort boxes in descending order
boxes.sort(reverse=True)

def can_lift(box, lifters, used):
    # Try to find a combination of lifters to lift the box
    for r in range(1, len(lifters) + 1):
        for combo in combinations(enumerate(lifters), r):
            indices, weights = zip(*combo)
            if sum(weights) >= box and all(not used[i] for i in indices):
                return indices
    return None

steps = []
while boxes and len(steps) < 6:
    step = []
    used_lifters = [False] * len(lifters)
    
    for box in boxes[:]:
        lifter_indices = can_lift(box, lifters, used_lifters)
        if lifter_indices is not None:
            step.append((box, list(lifter_indices)))
            for index in lifter_indices:
                used_lifters[index] = True
            boxes.remove(box)
    
    steps.append(step)

# Check if all boxes are lifted
if boxes:
    print("Unable to lift all boxes within 6 steps.")
else:
    # Print the steps
    output = []
    for i, step in enumerate(steps):
        output.append(f"Step {i+1}: {step}")
    print("<<<" + "\n".join(output) + ">>>")