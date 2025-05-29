from itertools import combinations

boxes = [356, 292, 240, 155, 142, 223, 178, 231, 254, 315, 236, 328, 203, 397, 44, 200, 80, 131, 175, 216]
lifters = [155, 116, 53, 98, 96, 156, 167]

# Sort boxes in descending order
boxes.sort(reverse=True)

steps = []
max_steps = 7

def find_lifter_combination(box, lifters, used):
    for r in range(1, len(lifters) + 1):
        for combo in combinations(range(len(lifters)), r):
            if all(not used[i] for i in combo):
                if sum(lifters[i] for i in combo) >= box:
                    return combo
    return None

for _ in range(max_steps):
    step = []
    used_lifters = [False] * len(lifters)
    remaining_boxes = []

    for box in boxes:
        lifter_indices = find_lifter_combination(box, lifters, used_lifters)
        if lifter_indices is not None:
            step.append((box, list(lifter_indices)))
            for index in lifter_indices:
                used_lifters[index] = True
        else:
            remaining_boxes.append(box)

    steps.append(step)
    boxes = remaining_boxes

    if not boxes:
        break

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")