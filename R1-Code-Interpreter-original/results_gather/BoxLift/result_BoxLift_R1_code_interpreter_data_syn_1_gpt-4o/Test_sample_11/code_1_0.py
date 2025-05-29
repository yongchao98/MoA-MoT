from itertools import combinations

boxes = [196, 90, 142, 22, 143, 101, 43, 166, 79, 129, 77, 187]
lifters = [71, 58, 104, 105]

# Sort boxes in descending order
boxes.sort(reverse=True)

steps = []

def find_lifter_combination(box_weight, lifters):
    # Try to find the smallest combination of lifters that can lift the box
    for r in range(1, len(lifters) + 1):
        for combo in combinations(lifters, r):
            if sum(combo) >= box_weight:
                return combo
    return None

while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    for i, box in enumerate(boxes):
        # Find a combination of lifters to lift the box
        available_lifters = [l for j, l in enumerate(lifters) if not used_lifters[j]]
        lifter_combo = find_lifter_combination(box, available_lifters)
        if lifter_combo:
            lifter_indices = [lifters.index(l) for l in lifter_combo]
            for idx in lifter_indices:
                used_lifters[idx] = True
            step.append((box, lifter_indices))
            boxes[i] = None  # Mark the box as lifted
    boxes = [b for b in boxes if b is not None]  # Remove lifted boxes
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")