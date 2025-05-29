from itertools import combinations

# Box weights and lifter capacities
boxes = [151, 393, 143, 179, 88, 173, 344, 212, 308, 258, 113, 261, 371, 335, 236, 258, 124, 262, 108, 209, 232, 128, 266, 232]
lifters = [84, 94, 87, 108, 198, 95, 121]

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
    available_lifters = lifters[:]
    for box in boxes[:]:
        combo = find_lifter_combination(box, available_lifters)
        if combo:
            step.append((box, [lifters.index(l) for l in combo]))
            boxes.remove(box)
            for l in combo:
                available_lifters.remove(l)
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")