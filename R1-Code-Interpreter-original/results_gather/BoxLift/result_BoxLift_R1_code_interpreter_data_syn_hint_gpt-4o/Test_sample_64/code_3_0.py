from itertools import combinations

boxes = [45, 139, 291, 189, 171, 265, 139, 312, 289, 394, 183, 229, 321, 250, 226, 335, 59, 127, 298, 339, 256, 161, 288, 179]
lifters = [52, 99, 86, 56, 92, 134]

# Sort boxes in descending order
boxes.sort(reverse=True)

def can_lift(box, lifters):
    for r in range(1, len(lifters) + 1):
        for combo in combinations(lifters, r):
            if sum(combo) >= box:
                return list(combo)
    return None

steps = []
while boxes:
    step = []
    remaining_lifters = lifters[:]
    for box in boxes[:]:
        lifter_combo = can_lift(box, remaining_lifters)
        if lifter_combo:
            step.append((box, lifter_combo))
            for lifter in lifter_combo:
                remaining_lifters.remove(lifter)
            boxes.remove(box)
    steps.append(step)
    if len(steps) > 14:
        break

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")