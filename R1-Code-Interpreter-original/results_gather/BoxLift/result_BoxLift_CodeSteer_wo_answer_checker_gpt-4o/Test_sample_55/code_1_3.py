boxes = [24, 78, 95, 23, 57, 36, 40, 26]
lifters = [56, 51, 66, 65]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []

def can_lift(box_weight, lifter_indices):
    return sum(lifters[i] for i in lifter_indices) >= box_weight

def find_lifters_for_box(box_weight, available_lifters):
    from itertools import combinations
    for r in range(1, len(available_lifters) + 1):
        for combo in combinations(available_lifters, r):
            if can_lift(box_weight, combo):
                return combo
    return []

def assign_boxes_to_lifters():
    remaining_boxes = boxes[:]
    for step in range(3):
        step_lift = []
        used_lifters = set()
        for box in remaining_boxes[:]:  # Use a copy of the list to modify it while iterating
            available_lifters = [i for i in range(len(lifters)) if i not in used_lifters]
            lifter_indices = find_lifters_for_box(box, available_lifters)
            if lifter_indices:
                step_lift.append((box, lifter_indices))
                used_lifters.update(lifter_indices)
                remaining_boxes.remove(box)
        steps.append(step_lift)
        if not remaining_boxes:
            break

assign_boxes_to_lifters()

# Format the output
output = ""
for i, step in enumerate(steps):
    output += f"Step {i + 1}: {[(box, list(lifter_indices)) for box, lifter_indices in step]}\n"

print(f"<<<{output.strip()}>>>")