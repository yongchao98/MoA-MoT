from itertools import combinations

boxes = [187, 344, 112, 66, 384, 247, 184, 298, 131, 51, 62, 255, 120, 357, 399, 287, 231, 161, 336, 256, 328, 239, 365, 245]
lifters = [165, 85, 52, 82, 114, 75]

# Sort boxes in descending order
boxes.sort(reverse=True)

def can_lift(box, lifters):
    for r in range(1, len(lifters) + 1):
        for combo in combinations(lifters, r):
            if sum(combo) >= box:
                return list(combo)
    return None

steps = []
while boxes and len(steps) < 11:
    step = []
    remaining_boxes = []
    used_lifters = [False] * len(lifters)
    
    for box in boxes:
        lifter_combo = can_lift(box, [l for i, l in enumerate(lifters) if not used_lifters[i]])
        if lifter_combo:
            step.append((box, lifter_combo))
            for lifter in lifter_combo:
                used_lifters[lifters.index(lifter)] = True
        else:
            remaining_boxes.append(box)
    
    if step:
        steps.append(step)
    boxes = remaining_boxes

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")