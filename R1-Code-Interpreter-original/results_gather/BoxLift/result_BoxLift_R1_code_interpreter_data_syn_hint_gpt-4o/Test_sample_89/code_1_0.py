from itertools import combinations

# Box weights and lifter capacities
boxes = [227, 106, 388, 173, 63, 178, 261, 205, 323, 124, 386, 81, 71, 127, 134, 212, 150, 114, 41, 277, 123, 152, 47, 47]
lifters = [166, 106, 122, 181, 102, 119, 147]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Function to find a combination of lifters to lift a box
def find_lifters_for_box(box_weight, available_lifters):
    for r in range(1, len(available_lifters) + 1):
        for combo in combinations(available_lifters, r):
            if sum(combo) >= box_weight:
                return combo
    return None

# Backtracking function to find the optimal steps
def backtrack(boxes, lifters, step, steps, max_steps):
    if not boxes:
        if len(steps) < max_steps[0]:
            max_steps[0] = len(steps)
            best_steps[:] = steps[:]
        return

    if len(steps) >= max_steps[0]:
        return

    available_lifters = lifters[:]
    current_step = []
    i = 0
    while i < len(boxes):
        box_weight = boxes[i]
        lifter_combo = find_lifters_for_box(box_weight, available_lifters)
        if lifter_combo:
            current_step.append((box_weight, [lifters.index(l) for l in lifter_combo]))
            for l in lifter_combo:
                available_lifters.remove(l)
            boxes.pop(i)
        else:
            i += 1

    steps.append(current_step)
    backtrack(boxes, lifters, step + 1, steps, max_steps)
    steps.pop()
    boxes.extend([b[0] for b in current_step])

# Initialize variables
best_steps = []
max_steps = [6]

# Start backtracking
backtrack(boxes, lifters, 0, [], max_steps)

# Print the best steps found
for idx, step in enumerate(best_steps):
    print(f"Step {idx + 1}: {step}")