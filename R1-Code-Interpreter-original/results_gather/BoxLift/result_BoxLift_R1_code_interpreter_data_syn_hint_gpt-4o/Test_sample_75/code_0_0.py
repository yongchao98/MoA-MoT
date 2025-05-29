from itertools import combinations

# Box weights and lifter capacities
boxes = [97, 83, 76, 46, 45, 36, 30, 14]
lifters = [69, 68, 59, 50]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

# Function to find the best combination of lifters for a given box
def find_lifters_for_box(box_weight, available_lifters):
    for r in range(1, len(available_lifters) + 1):
        for combo in combinations(available_lifters, r):
            if sum(combo) >= box_weight:
                return list(combo)
    return []

# Assign boxes to lifters in steps
steps = []
for _ in range(3):  # Maximum 3 steps
    step = []
    remaining_lifters = lifters[:]
    for box in boxes[:]:
        lifter_combo = find_lifters_for_box(box, remaining_lifters)
        if lifter_combo:
            step.append((box, lifter_combo))
            boxes.remove(box)
            for lifter in lifter_combo:
                remaining_lifters.remove(lifter)
    if step:
        steps.append(step)
    if not boxes:
        break

# Print the steps
for i, step in enumerate(steps, 1):
    print(f"Step {i}: {step}")