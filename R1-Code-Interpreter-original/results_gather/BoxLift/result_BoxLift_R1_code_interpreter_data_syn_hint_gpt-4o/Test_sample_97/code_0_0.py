from itertools import combinations

# Box weights and lifter capacities
boxes = [184, 93, 275, 216, 137, 181, 31, 79, 56, 138, 81, 205, 108, 193, 230, 252]
lifters = [47, 157, 156, 45, 151, 51]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

# Function to find the best combination of lifters to lift a box
def find_lifter_combination(box_weight, available_lifters):
    for r in range(1, len(available_lifters) + 1):
        for combo in combinations(available_lifters, r):
            if sum(combo) >= box_weight:
                return combo
    return None

# Assign boxes to lifters in steps
steps = []
while boxes:
    step = []
    used_lifters = set()
    for box in boxes[:]:
        available_lifters = [l for i, l in enumerate(lifters) if i not in used_lifters]
        combo = find_lifter_combination(box, available_lifters)
        if combo:
            step.append((box, [lifters.index(l) for l in combo]))
            used_lifters.update(lifters.index(l) for l in combo)
            boxes.remove(box)
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")