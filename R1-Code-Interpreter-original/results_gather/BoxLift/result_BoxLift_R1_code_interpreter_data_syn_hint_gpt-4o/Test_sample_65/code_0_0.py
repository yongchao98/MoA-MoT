from itertools import combinations

# Box weights and lifter capacities
boxes = [70, 51, 67, 51, 34, 82, 20, 93]
lifters = [64, 51, 79, 80]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

# Function to find the best combination of lifters for a given box
def find_lifter_combination(box_weight, available_lifters):
    for r in range(1, len(available_lifters) + 1):
        for combo in combinations(available_lifters, r):
            if sum(combo) >= box_weight:
                return combo
    return None

# Assign boxes to lifters in steps
steps = []
for _ in range(3):  # Maximum of 3 steps
    step = []
    available_lifters = lifters[:]
    for box in boxes[:]:
        combo = find_lifter_combination(box, available_lifters)
        if combo:
            step.append((box, [lifters.index(l) for l in combo]))
            boxes.remove(box)
            for l in combo:
                available_lifters.remove(l)
    if step:
        steps.append(step)
    if not boxes:
        break

# Print the steps
for i, step in enumerate(steps, 1):
    print(f"Step {i}: {step}")